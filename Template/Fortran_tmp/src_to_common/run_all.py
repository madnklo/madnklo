import argparse
import glob
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed



def check_onshell_widths():
    """Warn if the user has set a non-zero width for a particle that must be
    stable and on-shell in this process (e.g. WT for the top quark).
    The list of affected PDGs is recorded in Cards/onshell_zero_widths.dat
    at generation time.  If the file is absent the check is silently skipped.
    """
    import sys
    base = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    sentinel = os.path.join(base, 'Cards', 'onshell_zero_widths.dat')
    if not os.path.isfile(sentinel):
        return
    with open(sentinel) as f:
        pdgs = [int(l.strip()) for l in f if l.strip()]
    if not pdgs:
        return

    param_card = os.path.join(base, 'Cards', 'param_card.dat')
    if not os.path.isfile(param_card):
        return

    widths = {}
    in_decay = False
    with open(param_card) as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith('#'):
                continue
            upper = stripped.upper()
            if upper.startswith('BLOCK') or upper.startswith('DECAY'):
                in_decay = upper.startswith('DECAY')
            if in_decay:
                parts = stripped.split()
                if len(parts) >= 3 and parts[0].upper() == 'DECAY':
                    try:
                        widths[int(parts[1])] = float(parts[2])
                    except ValueError:
                        pass

    bad = [(pdg, widths[pdg]) for pdg in pdgs if pdg in widths and widths[pdg] != 0.0]
    if bad:
        use_color = hasattr(sys.stderr, 'isatty') and sys.stderr.isatty()
        RED  = '\033[1;31m' if use_color else ''
        BOLD = '\033[1m'    if use_color else ''
        RST  = '\033[0m'    if use_color else ''
        sep = RED + '=' * 70 + RST
        sys.stderr.write('\n' + sep + '\n')
        for pdg, val in bad:
            sys.stderr.write(
                RED + BOLD +
                'WARNING: Width of particle PDG=%d (WT) = %.6g in param_card.dat' % (pdg, val) +
                RST + '\n')
        sys.stderr.write(BOLD + 'This particle is a stable on-shell final state in this process.' + RST + '\n')
        sys.stderr.write(BOLD + 'The Fortran matrix elements assume WT=0; a non-zero value gives wrong results.' + RST + '\n')
        sys.stderr.write(BOLD + 'Please set WT=0 in Cards/param_card.dat before rerunning.' + RST + '\n')
        sys.stderr.write(sep + '\n\n')

def run(cmd, cwd):
    subprocess.run(cmd, shell=True, check=True, cwd=cwd)


# The 'libs' target of every P1_* makefile runs 'make' in the shared Source
# directories (the contribution's Source/DHELAS and the top-level Source). The
# top-level Source rebuilds lib/libpdf.a at every call (its 'all' target ends with
# 'rm -f PDF/*.o', which makes the PDF directory newer than the library), so
# 'make' run concurrently in several P1_* directories races on libpdf.a (e.g. an
# archive member deleted or truncated while another process is archiving or
# linking). Hence: build 'libs' once per P1_* directory, serially, before the
# parallel phase, which then builds only the executables.
def make_targets_without_libs(p):
    """Targets of the 'all' rule of p/makefile, except 'libs' (None if not found)."""
    try:
        with open(os.path.join(p, "makefile")) as f:
            text = f.read().replace("\\\n", " ")
    except IOError:
        return None
    for line in text.splitlines():
        if line.startswith("all:"):
            return [t for t in line[len("all:"):].split() if t != "libs"]
    return None


def make_libs_serially(paths):
    for p in paths:
        run("make libs", p)
        remove_outdated_executables(p)


# ---------------- cards -> compiled settings ----------------
# The run_card values are compiled in (Source/run_card.inc, included by
# setrun.f in lib/libgeneric.a), while param_card.dat and MadLoopParams.dat are
# read at run time. When Cards/run_card.dat is newer than Source/run_card.inc,
# run_card.inc is regenerated here (as 'bin/treatcards run' does, with the
# Python 2 MadGraph code); 'make libs' then rebuilds libgeneric.a, and the
# executables older than run_card.inc are removed so that make relinks them.
RUN_CARD_INC_SCRIPT = """
import sys
sys.path.insert(0, sys.argv[1])
import madgraph.various.banner as banner_mod
run_card = banner_mod.RunCardNLO(sys.argv[2])
if run_card['lpp1'] + run_card['lpp2'] == 1:
    run_card['lpp1'] = 0
    run_card['lpp2'] = 0
    run_card['ebeam1'] = 0
    run_card['ebeam2'] = 0
run_card.write_include_file(sys.argv[3])
"""


def output_dir():
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def find_python2():
    for exe in ("python2", "python2.7", "python"):
        path = shutil.which(exe)
        if not path:
            continue
        try:
            out = subprocess.run([path, "-c", "import sys; print(sys.version_info[0])"],
                                 capture_output=True, text=True, timeout=60)
        except Exception:
            continue
        if out.returncode == 0 and out.stdout.strip() == "2":
            return path
    return None


def run_card_include():
    return os.path.join(output_dir(), "Source", "run_card.inc")


def update_run_card_include():
    """Regenerate Source/run_card.inc if Cards/run_card.dat is newer."""
    card = os.path.join(output_dir(), "Cards", "run_card.dat")
    inc = run_card_include()
    if not os.path.isfile(card):
        return
    if os.path.isfile(inc) and os.path.getmtime(inc) >= os.path.getmtime(card):
        return
    python2 = find_python2()
    if python2 is None:
        sys.exit("Cards/run_card.dat was modified, but no Python 2 interpreter "
                 "(python2, python2.7 or python) was found to regenerate "
                 "Source/run_card.inc. Run 'bin/treatcards run' by hand.")
    print("Cards/run_card.dat was modified: regenerating Source/run_card.inc")
    subprocess.run([python2, "-c", RUN_CARD_INC_SCRIPT, os.path.dirname(output_dir()),
                    card, os.path.join(output_dir(), "Source")], check=True)
    if not (os.path.isfile(inc) and os.path.getmtime(inc) >= os.path.getmtime(card)):
        sys.exit("Source/run_card.inc was not regenerated: run 'bin/treatcards run' by hand.")


# Real emission: the sources, the executable and the outputs of each sector are
# in the sub-folder sector_<i>_<j> of the P1_* directory (the makefile of the
# P1_* directory builds sector_<i>_<j>/sector_<i>_<j>), and the executable is run
# there. Directories generated before this layout have everything in P1_*.
def exe_dir(p, exe):
    """Directory where executable exe of the P1_* directory p is and runs."""
    sub = os.path.join(p, exe)
    return sub if os.path.isdir(sub) else p


def remove_outdated_executables(p):
    """Remove the executables of p older than Source/run_card.inc, so that make
    relinks them with the rebuilt libgeneric.a (the makefiles do not list the
    libraries among the prerequisites of the executables)."""
    inc = run_card_include()
    if not os.path.isfile(inc):
        return
    t_inc = os.path.getmtime(inc)
    for target in make_targets_without_libs(p) or []:
        exe = os.path.join(exe_dir(p, target), target)
        if os.path.isfile(exe) and os.path.getmtime(exe) < t_inc:
            os.remove(exe)


def make_executables(p):
    """Build the executables of the 'all' rule of p (the libraries are built
    beforehand, serially, by make_libs_serially)."""
    targets = make_targets_without_libs(p)
    if targets:
        run("make " + " ".join(targets), p)
    else:
        run("make", p)
    return p


def collect_p1_paths(pattern):
    paths = []
    for path in sorted(glob.glob(pattern)):
        paths.extend(sorted(glob.glob(f"{path}/P1_*")))
    return paths


# ---------------- results bookkeeping ----------------
# Each executable writes results_<tag>.log and seed_<iseed>_plot_<tag>.dat in
# the directory where it runs (P1_*, or P1_*/sector_<i>_<j> for the real
# emission). The files produced by the last *complete* run of each
# contribution (B, V, R) are recorded in Events/run_manifest.json, and
# collect_results.py / collect_plots.py only use those: other files left in the
# P1_* directories (other seeds, saved copies, sectors that no longer exist)
# are ignored. The entry of a contribution is removed when its run starts and
# written back only if all its executables succeed.
def manifest_path():
    return os.path.join(output_dir(), "Events", "run_manifest.json")


def load_manifest():
    try:
        with open(manifest_path()) as f:
            manifest = json.load(f)
    except (IOError, ValueError):
        manifest = {}
    manifest.setdefault("format", 1)
    manifest.setdefault("contributions", {})
    return manifest


def save_manifest(manifest):
    path = manifest_path()
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(manifest, f, indent=1, sort_keys=True)
    os.replace(tmp, path)


def output_tag(exe):
    """'Born' -> 'B', 'virtual' -> 'V', 'sector_3_4' -> 'R_3_4'"""
    if exe == "Born":
        return "B"
    if exe == "virtual":
        return "V"
    if exe.startswith("sector_"):
        return "R_" + exe[len("sector_"):]
    raise ValueError(f"unknown executable '{exe}'")


def relpath(path):
    return os.path.relpath(os.path.abspath(path), output_dir())


def run_executable(task):
    """Run one executable and return the files it has written."""
    p, exe = task
    p = exe_dir(p, exe)
    start = time.time() - 1.0  # tolerance for the file-system time resolution
    run(f"./{exe}", p)
    tag = output_tag(exe)
    results = os.path.join(p, f"results_{tag}.log")
    if not (os.path.isfile(results) and os.path.getmtime(results) >= start):
        raise RuntimeError(f"{exe} in {p} did not write results_{tag}.log")
    plots = [f for f in glob.glob(os.path.join(p, f"seed_*_plot_{tag}.dat"))
             if os.path.getmtime(f) >= start]
    if not plots:
        raise RuntimeError(f"{exe} in {p} did not write seed_*_plot_{tag}.dat")
    plot = max(plots, key=os.path.getmtime)
    return {"dir": relpath(p), "executable": exe,
            "results": relpath(results), "results_mtime": os.path.getmtime(results),
            "plot": relpath(plot), "plot_mtime": os.path.getmtime(plot)}


def pool_map(fn, items, workers, label):
    """Apply fn to items with at most 'workers' processes; wait for all of them,
    then raise if any failed. Returns the results in the order of items."""
    results = [None] * len(items)
    failures = []
    with ProcessPoolExecutor(max_workers=max(1, workers)) as executor:
        futures = {executor.submit(fn, item): i for i, item in enumerate(items)}
        for future in as_completed(futures):
            i = futures[future]
            name = items[i] if isinstance(items[i], str) else os.path.join(*items[i])
            try:
                results[i] = future.result()
                print(f"[{label} OK] {name}")
            except Exception as e:
                failures.append(name)
                print(f"[{label} FAIL] {name}: {e}")
    if failures:
        raise RuntimeError(f"{len(failures)} {label} task(s) failed: " + ", ".join(failures))
    return results


def run_contribution(label, pattern, n, compile_jobs):
    """Libraries (serially), executables (compile_jobs P1_* directories at a
    time), then all the executables of the contribution with a single pool of
    n processes (e.g. all the sectors of all the real-emission directories)."""
    paths = collect_p1_paths(pattern)
    if not paths:
        print(f"No P1_* directories found for {pattern}.")
        return
    manifest = load_manifest()
    manifest["contributions"].pop(label, None)
    save_manifest(manifest)

    print(f"--- {label}: libraries ({len(paths)} P1_* directories, serially)")
    make_libs_serially(paths)
    print(f"--- {label}: compiling, {compile_jobs} P1_* directories at a time")
    pool_map(make_executables, paths, compile_jobs, "compile")

    tasks = [(p, exe) for p in paths for exe in (make_targets_without_libs(p) or [])]
    print(f"--- {label}: running {len(tasks)} executables, {n} at a time")
    runs = pool_map(run_executable, tasks, n, "run")

    manifest = load_manifest()
    manifest["contributions"][label] = {
        "finished": time.strftime("%Y-%m-%d %H:%M:%S"),
        "runs": runs,
    }
    save_manifest(manifest)
    print(f"--- {label}: done, {len(runs)} results recorded in {relpath(manifest_path())}")


# ---------------- contributions ----------------
CONTRIBUTIONS = {
    "B": "../LO_*/SubProcesses",
    "V": "../NLO_V_*/SubProcesses",
    "R": "../NLO_R_*/SubProcesses",
}
MODES = {
    "B": ["B"],
    "V": ["V"],
    "R": ["R"],
    "NLOcorr": ["V", "R"],
    "ALL": ["B", "V", "R"],
}
NAMES = {"B": "LO (Born)", "V": "NLO_V (Virtual)", "R": "NLO_R (Real)"}


# ---------------- test-only mode ----------------
# 'run_all.py <mode> --test_only' compiles as usual, then runs each virtual and
# real-emission executable with '--test_only' (which can also be done by hand in
# the P1_* directory, './virtual --test_only', or in the sector folder,
# 'cd sector_3_4; ./sector_3_4 --test_only'):
# the executables only write the pole test of V+I (test_poles_V.log) and the
# IR-limit tests of R-K (testR_<i>_<j>.log), then stop, without integrating and
# without touching results_*.log, the plots or Events/run_manifest.json.
def test_log_name(exe):
    if exe == "virtual":
        return "test_poles_V.log"
    if exe.startswith("sector_"):
        return f"testR_{exe[len('sector_'):]}.log"
    return None


def test_executable(task):
    p, exe = task
    p = exe_dir(p, exe)
    start = time.time() - 1.0
    run(f"./{exe} --test_only", p)
    log = os.path.join(p, test_log_name(exe))
    if not (os.path.isfile(log) and os.path.getmtime(log) >= start):
        raise RuntimeError(f"{exe} in {p} did not write {test_log_name(exe)}")
    return log


def summarize_poles(log):
    """Relative |V+I|/max(|V|,|I|) of the double and single poles, per point."""
    rel = {"Double": [], "Single": []}
    pattern = re.compile(r"\s*(Double|Single) pole V, I, sum\s+(\S+)\s+(\S+)\s+(\S+)")
    with open(log) as f:
        for line in f:
            m = pattern.match(line)
            if not m:
                continue
            v, i, s = (float(x.replace("D", "E")) for x in m.group(2, 3, 4))
            scale = max(abs(v), abs(i))
            rel[m.group(1)].append(abs(s) / scale if scale > 0 else abs(s))
    return rel


# Limits tested by testR_<i>_<j>.f: Si, Sj, Cij are singular limits of the
# sector, where R and the counterterm K must agree (|R-K|/|K| -> 0); Cir and
# Cjr are spurious limits (damped by the sector function), where R and K must
# only stay finite.
SPURIOUS_LIMITS = ("Cir", "Cjr")


def summarize_limits(log):
    """{limit: [rows of one test event, ...]}, rows = (lambda, R, K)."""
    res = {}
    lim, rows = None, []

    def flush():
        if lim is not None and rows:
            res.setdefault(lim, []).append(rows)

    with open(log) as f:
        for line in f:
            if "EVENT NUMBER" in line or line.strip().startswith("LIM ="):
                flush()
                lim = line.split("=", 1)[1].strip() if line.strip().startswith("LIM =") else None
                rows = []
                continue
            parts = line.split()
            if lim is None or len(parts) < 3:
                continue
            try:
                lam, r, k = (float(x.replace("D", "E")) for x in parts[:3])
            except ValueError:
                continue
            rows.append((lam, r, k))
    flush()
    return res


def check_limit(lim, events, tol):
    """(value, ok) for one limit over all its test events. Singular limits: the
    best |R-K|/|K| along the approach, worst event. Spurious limits: the largest
    growth of |R| or |K| from lambda ~ 1e-2 to the smallest lambda (finite: ~1)."""
    spurious = lim.startswith(SPURIOUS_LIMITS)
    worst = 0.0
    for rows in events:
        if spurious:
            ref = min(rows, key=lambda row: abs(math.log10(abs(row[0]) or 1e-300) + 2))
            last = rows[-1]
            growth = max(abs(last[1]) / abs(ref[1]) if ref[1] else 0.0,
                         abs(last[2]) / abs(ref[2]) if ref[2] else 0.0)
            worst = max(worst, growth)
        else:
            rels = [abs(r - k) / abs(k) for _, r, k in rows if k != 0.0]
            worst = max(worst, min(rels) if rels else float("inf"))
    if spurious:
        return worst, worst <= 10.0
    return worst, worst <= tol


def test_header(title, tol_poles, tol_limits):
    return [f"{title}, {time.strftime('%Y-%m-%d %H:%M:%S')}",
            f"tolerances: poles {tol_poles:g} (|V+I|/max(|V|,|I|)), "
            f"limits {tol_limits:g} (best |R-K|/|K| along the approach, worst event; "
            f"spurious limits Cir, Cjr: R and K must stay finite)", ""]


def summarize_test_log(exe, log, tol_poles, tol_limits, show_time=False):
    """(line, ok) summarizing the test log of executable exe ('virtual' or
    'sector_<i>_<j>')."""
    where = relpath(log)
    if show_time:
        where += " (" + time.strftime("%Y-%m-%d %H:%M", time.localtime(os.path.getmtime(log))) + ")"
    if exe == "virtual":
        rel = summarize_poles(log)
        if not rel["Double"]:
            return f"[CHECK] {where}: no pole tests found", False
        worst = max(rel["Double"] + rel["Single"])
        ok = worst <= tol_poles
        return (f"[{'OK' if ok else 'CHECK':5}] {where}: {len(rel['Double'])} points, max relative "
                f"V+I: double pole {max(rel['Double']):.2e}, "
                f"single pole {max(rel['Single'] + [0.0]):.2e}"), ok
    res = summarize_limits(log)
    if not res:
        return f"[CHECK] {where}: no limit tests found", False
    parts = []
    all_ok = True
    for lim, events in res.items():
        value, ok = check_limit(lim, events, tol_limits)
        all_ok = all_ok and ok
        if lim.startswith(SPURIOUS_LIMITS):
            parts.append(f"{lim} " + ("finite" if ok else f"grows x{value:.1e}(!)"))
        else:
            parts.append(f"{lim} {value:.1e}" + ("" if ok else "(!)"))
    return f"[{'OK' if all_ok else 'CHECK':5}] {where}: " + ", ".join(parts), all_ok


def limit_details(log):
    """Approach to the limits of one testR log: one table per limit, with,
    for each lambda, the worst and the median of |R-K|/|K| over the test
    events (8 significant digits)."""
    res = summarize_limits(log)
    if not res:
        return []
    lines = [relpath(log)]
    for lim, events in res.items():
        table = {}
        for rows in events:
            for lam, r, k in rows:
                if not lam:
                    continue
                key = round(math.log10(abs(lam)), 3)
                rel = abs(r - k) / abs(k) if k != 0.0 else float("inf")
                table.setdefault(key, []).append(rel)
        kind = ("spurious limit: R and K are only required to stay finite"
                if lim.startswith(SPURIOUS_LIMITS) else "singular limit: |R-K|/|K| -> 0")
        lines.append(f"  limit {lim} ({kind})")
        lines.append(f"    |R-K|/|K| at each lambda, over the {len(events)} test events")
        lines.append(f"    {'lambda':>8}   {'worst':>15}   {'median':>15}")
        lines.append("    " + "-" * 44)
        for key in sorted(table, reverse=True):
            vals = sorted(table[key])
            lines.append(f"    {10 ** key:8.0e}   {vals[-1]:15.7e}   {vals[len(vals) // 2]:15.7e}")
        lines.append("")
    return lines


def details_section(entries):
    """Detailed approach to the limits for the real-emission test logs
    entries = [(exe, log), ...]; appended after the summary."""
    lines = []
    for exe, log in entries:
        if exe.startswith("sector_"):
            lines.extend(limit_details(log))
    if lines:
        lines = ["", "=" * 78, "Details: approach to the limits, lambda -> 0", "=" * 78, ""] + lines
    return lines


def run_tests(labels, n, compile_jobs, tol_poles, tol_limits):
    lines = test_header("run_all.py --test_only", tol_poles, tol_limits)
    n_check = 0
    entries = []
    for label in labels:
        if label == "B":
            print("--- B: nothing to test")
            continue
        paths = collect_p1_paths(CONTRIBUTIONS[label])
        if not paths:
            continue
        print(f"--- {label}: libraries ({len(paths)} P1_* directories, serially)")
        make_libs_serially(paths)
        print(f"--- {label}: compiling, {compile_jobs} P1_* directories at a time")
        pool_map(make_executables, paths, compile_jobs, "compile")
        tasks = [(p, exe) for p in paths for exe in (make_targets_without_libs(p) or [])
                 if test_log_name(exe)]
        print(f"--- {label}: testing {len(tasks)} executables, {n} at a time")
        logs = pool_map(test_executable, tasks, n, "test")
        for (p, exe), log in zip(tasks, logs):
            line, ok = summarize_test_log(exe, log, tol_poles, tol_limits)
            n_check += not ok
            lines.append(line)
            entries.append((exe, log))
        lines.append("")
    lines.append(f"{n_check} test log(s) to check" if n_check else "all tests within tolerance")
    details = details_section(entries)

    summary = next_free_name(os.path.join(output_dir(), "Events", "tests"), ".txt")
    os.makedirs(os.path.dirname(summary), exist_ok=True)
    with open(summary, "w") as f:
        f.write("\n".join(lines + details) + "\n")
    print("\n" + "\n".join(lines))
    print(f"\nTest summary (with the approach to the limits) written to {relpath(summary)}")


def collect_test_logs(tol_poles, tol_limits):
    """Summary of the test logs currently present (from 'run_all.py --test_only',
    from the executables run by hand with --test_only, or from the first points
    of a normal run): one line per virtual and per sector of the current
    makefiles, with the time each log was written; missing logs are listed."""
    lines = test_header("summary of the test logs", tol_poles, tol_limits)
    n_check = 0
    entries = []
    for label in ("V", "R"):
        missing = []
        for p in collect_p1_paths(os.path.join(output_dir(), CONTRIBUTIONS[label][3:])):
            for exe in make_targets_without_libs(p) or []:
                name = test_log_name(exe)
                if not name:
                    continue
                log = os.path.join(exe_dir(p, exe), name)
                if not os.path.isfile(log):
                    missing.append(relpath(log))
                    continue
                line, ok = summarize_test_log(exe, log, tol_poles, tol_limits, show_time=True)
                n_check += not ok
                lines.append(line)
                entries.append((exe, log))
        for m in missing:
            lines.append(f"[MISSING] {m}")
        n_check += len(missing)
        lines.append("")
    lines.append(f"{n_check} test log(s) to check" if n_check else "all tests within tolerance")
    return lines + details_section(entries)


def next_free_name(stem, ext=""):
    i = 1
    while os.path.exists(f"{stem}_{i}{ext}"):
        i += 1
    return f"{stem}_{i}{ext}"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run LO/NLO computations in parallel across P1_* directories"
    )
    parser.add_argument(
        "mode",
        choices=list(MODES),
        help="B = Born, V = Virtual, R = Real, NLOcorr = V+R, ALL = B+V+R",
    )
    parser.add_argument(
        "-n",
        type=int,
        default=os.cpu_count(),
        help="Number of executables (Born, virtual, sector_*) run at the same time",
    )
    parser.add_argument(
        "--test_only",
        action="store_true",
        help="Only run the tests: pole cancellation in the virtual, IR limits in the "
             "real-emission sectors (summary in Events/tests_N.txt); no integration",
    )
    parser.add_argument(
        "--tol-poles",
        type=float,
        default=1e-5,
        help="--test_only: tolerance on |V+I|/max(|V|,|I|) for the poles (default 1e-5)",
    )
    parser.add_argument(
        "--tol-limits",
        type=float,
        default=1e-3,
        help="--test_only: tolerance on the best |R-K|/|K| reached in each limit "
             "(default 1e-3)",
    )
    parser.add_argument(
        "-j", "--compile-jobs",
        type=int,
        default=2,
        help="Number of P1_* directories compiled at the same time (default 2; "
             "compilation, especially of the virtuals, can take a lot of memory)",
    )

    check_onshell_widths()
    args = parser.parse_args()
    # run from bin/, where the patterns above are relative to
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    update_run_card_include()

    try:
        if args.test_only:
            run_tests(MODES[args.mode], args.n, args.compile_jobs, args.tol_poles,
                      args.tol_limits)
            sys.exit(0)
        for label in MODES[args.mode]:
            print(f"=== Running {NAMES[label]} ===")
            run_contribution(label, CONTRIBUTIONS[label], args.n, args.compile_jobs)
    except (RuntimeError, subprocess.CalledProcessError) as e:
        sys.exit(f"\nrun_all.py stopped: {e}")
