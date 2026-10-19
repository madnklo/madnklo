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
import heapq
import itertools
from collections import deque
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED



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

# ---------------- output of make and of the executables ----------------
# By default the output of make goes to P1_*/compile.log and the output of each
# executable to log_<executable>.txt (log_<executable>_test.txt with
# --test_only) in the directory where it runs, and run_all.py prints one line
# per directory compiled and per executable run; with --verbose everything is
# printed on the screen, as before. The setting is passed to the worker
# processes through the environment (macOS starts them with 'spawn').
VERBOSE_ENV = "MADNKLO_RUN_ALL_VERBOSE"


def verbose():
    return bool(os.environ.get(VERBOSE_ENV))


ERROR_RE = re.compile(r"error|fatal|\*\*\*|undefined|not found|cannot|no such", re.I)


def tail(path, n=12, width=160):
    """The interesting end of a log, for an error message: the lines of the last
    command (after the last '###' header) around those that look like errors
    (or else the last n lines), each cut to width characters."""
    try:
        with open(path, errors="replace") as f:
            lines = f.read().splitlines()
    except IOError:
        return ""
    for i in range(len(lines) - 1, -1, -1):
        if lines[i].startswith("### "):
            lines = lines[i + 1:]
            break
    hits = [i for i, line in enumerate(lines) if ERROR_RE.search(line)]
    if hits:
        keep = sorted({j for i in hits for j in range(max(0, i - 3), min(len(lines), i + 2))})
        shown = [lines[j] for j in keep][-n:]
    else:
        shown = lines[-n:]
    return "\n".join("      | " + (line if len(line) <= width else line[:width] + " ...")
                     for line in shown)


def run(cmd, cwd, log=None, append=False):
    """Run cmd in cwd, its output in the file log (on the screen if log is None
    or with --verbose). On failure the last lines of log are in the error."""
    if log is None or verbose():
        subprocess.run(cmd, shell=True, check=True, cwd=cwd)
        return
    with open(log, "a" if append else "w") as f:
        f.write(f"### {time.strftime('%Y-%m-%d %H:%M:%S')}  {cmd}  (in {relpath(cwd)})\n")
        f.flush()
        returncode = subprocess.run(cmd, shell=True, cwd=cwd, stdout=f,
                                    stderr=subprocess.STDOUT).returncode
    if returncode != 0:
        raise RuntimeError(f"'{cmd}' failed (exit status {returncode}); "
                           f"last lines of {relpath(log)}:\n{tail(log)}")


def short(path):
    """How a directory is shown: relative to the output directory, without
    the 'SubProcesses/' level."""
    return relpath(path).replace(os.sep + "SubProcesses" + os.sep, os.sep)


def compile_log(p):
    return os.path.join(p, "compile.log")


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
    width = max(len(short(p)) for p in paths)
    for p in paths:
        if not verbose():
            print(f"    libraries  {short(p):<{width}}  ", end="", flush=True)
        start = time.time()
        try:
            run("make libs", p, compile_log(p))
        except RuntimeError:
            if not verbose():
                print("FAILED")
            raise
        if not verbose():
            print(f"{time.time() - start:7.1f} s")
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
        run("make " + " ".join(targets), p, compile_log(p), append=True)
    else:
        run("make", p, compile_log(p), append=True)
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
# Events/collect.py only uses those: other files left in the
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
    run(f"./{exe}", p, os.path.join(p, f"log_{exe}.txt"))
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


def timed_call(fn, item):
    start = time.time()
    return fn(item), time.time() - start


# ---------------- scheduling ----------------
# All the contributions of the mode share one set of n process slots. The
# libraries are built first, serially, for all the P1_* directories (they
# rebuild the shared Source, see make_libs_serially). Then the P1_* directories
# are compiled, at most compile_jobs at the same time, and each executable
# starts as soon as its directory is compiled and a slot is free: free slots go
# to compilation first, then to the executables in the order V, R, B (the
# virtuals are usually the slowest). So, e.g., the real-emission sectors run
# while the virtuals are still compiling or running. A failure in one
# contribution does not stop the others; each contribution is recorded in the
# manifest as soon as all its executables have succeeded.
# --sequential runs the contributions one after the other, as before.
RUN_ORDER = ("V", "R", "B")


def schedule(labels, n, compile_jobs, run_fn, run_verb, exe_filter=None, label_done=None):
    """Build and run the executables of the contributions labels (see above).
    run_fn((p, exe)) runs one executable; exe_filter(exe) selects the
    executables (default all); label_done(label, state) is called when all
    the tasks of a contribution are over. Returns {label: state}, state =
    {"tasks": [(p, exe), ...], "results": {(p, exe): result}, "failed": [names]}."""
    order = sorted(labels, key=RUN_ORDER.index)
    paths = {}
    for label in order:
        paths[label] = collect_p1_paths(CONTRIBUTIONS[label])
        if not paths[label]:
            print(f"--- {label}: no P1_* directories found for {CONTRIBUTIONS[label]}")
    order = [label for label in order if paths[label]]
    all_paths = [p for label in order for p in paths[label]]
    if not all_paths:
        return {}
    where = "" if verbose() else " (make output in P1_*/compile.log)"
    print(f"--- libraries, {len(all_paths)} P1_* directories one at a time{where}")
    make_libs_serially(all_paths)

    state = {}
    for label in order:
        tasks = [(p, exe) for p in paths[label] for exe in (make_targets_without_libs(p) or [])
                 if exe_filter is None or exe_filter(exe)]
        state[label] = {"tasks": tasks, "results": {}, "failed": [],
                        "left": len(paths[label]) + len(tasks)}
    n_tasks = sum(len(state[label]["tasks"]) for label in order)
    slots, compile_slots = max(1, n), max(1, min(n, compile_jobs))
    output = "" if verbose() else f", output in log_<executable>{'_test' if run_verb == 'tested' else ''}.txt"
    print(f"--- compiling {len(all_paths)} P1_* directories ({compile_slots} at a time) and "
          f"running their {n_tasks} executables as soon as they are compiled: "
          f"{slots} processes at a time, order {', '.join(order)}{output}")

    names = [short(p) for p in all_paths] + [short(os.path.join(*t)) for label in order
                                             for t in state[label]["tasks"]]
    width = max(len(name) for name in names)
    counters = {"compile": [0, len(all_paths)], "run": [0, n_tasks]}
    compile_queue = deque((label, p) for label in order for p in paths[label])
    ready, seq = [], itertools.count()
    running, n_compiling = {}, 0

    def finished(label):
        st = state[label]
        st["left"] -= 1
        if st["left"] == 0 and label_done is not None:
            label_done(label, st)

    with ProcessPoolExecutor(max_workers=slots) as executor:
        while compile_queue or ready or running:
            while len(running) < slots:
                if compile_queue and n_compiling < compile_slots:
                    label, p = compile_queue.popleft()
                    future = executor.submit(timed_call, make_executables, p)
                    running[future] = ("compile", label, p)
                    n_compiling += 1
                elif ready:
                    _, _, label, task = heapq.heappop(ready)
                    future = executor.submit(timed_call, run_fn, task)
                    running[future] = ("run", label, task)
                else:
                    break
            done, _ = wait(running, return_when=FIRST_COMPLETED)
            for future in done:
                kind, label, item = running.pop(future)
                st = state[label]
                if kind == "compile":
                    n_compiling -= 1
                    name, verb = short(item), "compiled"
                else:
                    name, verb = short(os.path.join(*item)), run_verb
                counter = counters[kind]
                counter[0] += 1
                k = f"{counter[0]:>{len(str(n_tasks if kind == 'run' else len(all_paths)))}}/{counter[1]}"
                try:
                    result, seconds = future.result()
                except Exception as e:
                    print(f"    [{k}] FAILED    {name}: {e}", flush=True)
                    st["failed"].append(name)
                    if kind == "compile":
                        # its executables cannot run
                        skipped = [t for t in st["tasks"] if t[0] == item]
                        counters["run"][1] -= len(skipped)
                        st["left"] -= len(skipped)
                    finished(label)
                    continue
                print(f"    [{k}] {verb:<9} {name:<{width}}  {seconds:7.1f} s", flush=True)
                if kind == "compile":
                    for task in st["tasks"]:
                        if task[0] == item:
                            heapq.heappush(ready, (RUN_ORDER.index(label), next(seq), label, task))
                else:
                    st["results"][item] = result
                finished(label)
    return state


def record_contribution(label, st):
    """Manifest entry of a contribution whose executables all succeeded."""
    if st["failed"]:
        print(f"--- {label}: {len(st['failed'])} failure(s), not recorded in "
              f"{relpath(manifest_path())}", flush=True)
        return
    manifest = load_manifest()
    manifest["contributions"][label] = {
        "finished": time.strftime("%Y-%m-%d %H:%M:%S"),
        "runs": [st["results"][task] for task in st["tasks"]],
    }
    save_manifest(manifest)
    print(f"--- {label}: done, {len(st['tasks'])} results recorded in "
          f"{relpath(manifest_path())}", flush=True)


def run_mode(labels, n, compile_jobs, sequential):
    """Run the contributions labels: all together (default) or one after the
    other (sequential, stopping at the first failure)."""
    manifest = load_manifest()
    for label in labels:
        manifest["contributions"].pop(label, None)
    save_manifest(manifest)
    groups = [[label] for label in labels] if sequential else [list(labels)]
    failed = []
    for group in groups:
        print(f"=== Running {' + '.join(NAMES[label] for label in group)}"
              + (" (in parallel)" if len(group) > 1 else "") + " ===")
        state = schedule(group, n, compile_jobs, run_executable, "done",
                         label_done=record_contribution)
        failed += [label for label in group if state.get(label, {}).get("failed")]
        if failed:
            break
    if failed:
        recorded = [label for label in labels if label in load_manifest()["contributions"]]
        raise RuntimeError(
            f"failures in {', '.join(NAMES[label] for label in failed)} (see the FAILED "
            f"lines above)" + (f"; recorded: {', '.join(recorded)} (collect them with "
                               f"'python3 Events/collect.py <MODE>')" if recorded else ""))


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
    run(f"./{exe} --test_only", p, os.path.join(p, f"log_{exe}_test.txt"))
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


def run_tests(mode, labels, n, compile_jobs, tol_poles, tol_limits):
    if not any(label in ("V", "R") for label in labels):
        print(f"--- {mode}: nothing to test (tests exist for V and R only)")
        return
    lines = test_header(f"run_all.py {mode} --test_only", tol_poles, tol_limits)
    n_check = 0
    entries = []
    test_labels = [label for label in labels if label in ("V", "R")]
    if "B" in labels:
        print("--- B: nothing to test")
    print(f"=== Testing {' + '.join(NAMES[label] for label in test_labels)} ===")
    state = schedule(test_labels, n, compile_jobs, test_executable, "tested",
                     exe_filter=test_log_name)
    failed = [name for st in state.values() for name in st["failed"]]
    if failed:
        raise RuntimeError(f"{len(failed)} failure(s): " + ", ".join(failed))
    for label in test_labels:
        if label not in state:
            continue
        for task in state[label]["tasks"]:
            p, exe = task
            log = state[label]["results"][task]
            line, ok = summarize_test_log(exe, log, tol_poles, tol_limits)
            n_check += not ok
            lines.append(line)
            entries.append((exe, log))
        lines.append("")
    lines.append(f"{n_check} test log(s) to check" if n_check else "all tests within tolerance")
    details = details_section(entries)

    summary = os.path.join(new_run_dir(f"run_{mode}_tests"), "tests.txt")
    with open(summary, "w") as f:
        f.write("\n".join(lines + details) + "\n")
    print("\n" + "\n".join(lines))
    print(f"\nTest summary (with the approach to the limits) written to {relpath(summary)}")


def collect_test_logs(tol_poles, tol_limits, labels=("V", "R")):
    """Summary of the test logs currently present (from 'run_all.py --test_only',
    from the executables run by hand with --test_only, or from the first points
    of a normal run): one line per virtual and per sector of the current
    makefiles, with the time each log was written; missing logs are listed."""
    lines = test_header("summary of the test logs", tol_poles, tol_limits)
    n_check = 0
    entries = []
    for label in [l for l in ("V", "R") if l in labels]:
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


def new_run_dir(stem):
    """Create and return Events/<stem>_<N>, with N = 1, 2, ... the first free
    number (as MG5's Events/run_NN; mkdir reserves the name, so that two runs
    started at the same time never share a directory)."""
    events = os.path.join(output_dir(), "Events")
    os.makedirs(events, exist_ok=True)
    i = 1
    while True:
        path = os.path.join(events, f"{stem}_{i}")
        try:
            os.mkdir(path)
            return path
        except FileExistsError:
            i += 1


# ---------------- collection ----------------
# At the end of a run, Events/collect.py <MODE> writes Events/run_<MODE>_<N>/
# (cross sections with scale variations, test summary, histograms and plots,
# copies of the cards); it can also be run by hand.
def collect_run(mode, tol_poles, tol_limits):
    script = os.path.join(output_dir(), "Events", "collect.py")
    if not os.path.isfile(script):
        print(f"\n{relpath(script)} not found: results not collected")
        return
    print("\n=== Collecting results ===")
    result = subprocess.run([sys.executable, script, mode,
                             f"--tol-poles={tol_poles}", f"--tol-limits={tol_limits}"])
    if result.returncode != 0:
        print(f"\nWARNING: {relpath(script)} {mode} failed (exit status "
              f"{result.returncode}); the results of the run are recorded in "
              f"{relpath(manifest_path())}: run it again by hand")


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
        help="Maximum number of processes at the same time, compilations and "
             "executables (Born, virtual, sector_*) of all the contributions of MODE "
             "together (default: the number of CPU cores, here %(default)s)",
    )
    parser.add_argument(
        "--test_only",
        action="store_true",
        help="Only run the tests: pole cancellation in the virtual, IR limits in the "
             "real-emission sectors (summary in Events/run_<MODE>_tests_N/tests.txt); "
             "no integration",
    )
    parser.add_argument(
        "--tol-poles",
        type=float,
        default=1e-5,
        help="Tests (--test_only, and tests.txt of the collection at the end of a run): "
             "tolerance on |V+I|/max(|V|,|I|) for the poles (default 1e-5)",
    )
    parser.add_argument(
        "--tol-limits",
        type=float,
        default=1e-3,
        help="Tests (--test_only, and tests.txt of the collection at the end of a run): "
             "tolerance on the best |R-K|/|K| reached in each limit (default 1e-3)",
    )
    parser.add_argument(
        "-j", "--compile-jobs",
        type=int,
        default=2,
        help="Of these, maximum number of P1_* directories compiled at the same time "
             "(default 2; compilation, especially of the virtuals, can take a lot of memory)",
    )
    parser.add_argument(
        "--sequential",
        action="store_true",
        help="Run the contributions of MODE one after the other (B, V, R), each "
             "compiled and run completely before the next one starts; by default "
             "they share the n processes, and an executable starts as soon as its "
             "directory is compiled (virtuals first)",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Print the full output of make and of the executables (by default it "
             "goes to P1_*/compile.log and log_<executable>.txt, and run_all.py "
             "prints one line per directory compiled and per executable run)",
    )
    parser.add_argument(
        "--no-collect",
        action="store_true",
        help="Do not run Events/collect.py at the end (by default it writes the "
             "results, tests and plots of the run to Events/run_<MODE>_N/)",
    )

    check_onshell_widths()
    args = parser.parse_args()
    if args.verbose:
        os.environ[VERBOSE_ENV] = "1"
    else:
        os.environ.pop(VERBOSE_ENV, None)
    # run from bin/, where the patterns above are relative to
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    update_run_card_include()

    try:
        if args.test_only:
            run_tests(args.mode, MODES[args.mode], args.n, args.compile_jobs,
                      args.tol_poles, args.tol_limits)
            sys.exit(0)
        run_mode(MODES[args.mode], args.n, args.compile_jobs, args.sequential)
    except (RuntimeError, subprocess.CalledProcessError) as e:
        sys.exit(f"\nrun_all.py stopped: {e}")
    if not args.no_collect:
        collect_run(args.mode, args.tol_poles, args.tol_limits)
