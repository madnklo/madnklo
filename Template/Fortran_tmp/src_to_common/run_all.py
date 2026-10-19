import argparse
import glob
import json
import os
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


def remove_outdated_executables(p):
    """Remove the executables of p older than Source/run_card.inc, so that make
    relinks them with the rebuilt libgeneric.a (the makefiles do not list the
    libraries among the prerequisites of the executables)."""
    inc = run_card_include()
    if not os.path.isfile(inc):
        return
    t_inc = os.path.getmtime(inc)
    for target in make_targets_without_libs(p) or []:
        exe = os.path.join(p, target)
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
# its P1_* directory. The files produced by the last *complete* run of each
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
        for label in MODES[args.mode]:
            print(f"=== Running {NAMES[label]} ===")
            run_contribution(label, CONTRIBUTIONS[label], args.n, args.compile_jobs)
    except (RuntimeError, subprocess.CalledProcessError) as e:
        sys.exit(f"\nrun_all.py stopped: {e}")
