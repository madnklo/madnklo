#!/usr/bin/env python3
"""Collect the outputs of bin/run_all.py into one directory per run.

    python3 Events/collect.py [MODE]        (MODE = B, V, R, NLOcorr, ALL; default ALL)

writes Events/run_<MODE>_<N>/ (N = 1, 2, ...: the first free number), with
  results.txt   cross sections of the contributions of MODE and of their sums
                (V + R, B + V + R), with the renormalisation-scale variations;
  tests.txt     summary of the test logs (V+I poles, IR limits of the sectors),
                for the V and R contributions of MODE;
  <X>.HwU       summed histograms, X = B, V, R, NLOcorr (V+R), NLO (B+V+R),
                with the scale-variation weights, and <X>.gnuplot / .pdf for
                the top combination(s) of MODE;
  sources.txt   the files used, and when the run_all.py runs finished;
  proc_card.dat, run_card.dat, param_card.dat, damping_factors.inc
                (and MadLoopParams.dat if V is collected)   copies of the cards:
                proc_card.dat has the MG5 commands that generated the
                output (model, process, subtraction scheme, options);
  <analysis>.f  copy of the analysis file (Common_Files/, 'analysis' in
                the run_card).
bin/run_all.py MODE runs it at the end of every run (unless --no-collect).
"""

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
from dataclasses import dataclass
from typing import Dict, List, Tuple


# ---------------------------------------------------------------------------
# cross sections in results_*.log
FLOAT = r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)"

SIGMA_RE = re.compile(
    r"sigma\s+\w+\s+\[pb\]\s*=\s*" + FLOAT + r"\s*\+\-\s*" + FLOAT
)

# scale variations (rw_rscale, rw_fscale in run_card.dat):
#   ' sigma <tag> [pb] muR/muR0 =   2.00000 :  value +- error'
# and, with factorisation-scale variations (hadron beams),
#   ' sigma <tag> [pb] muR/muR0 =   2.00000 muF/muF0 =   0.50000 :  value +- error'
SIGMA_VAR_RE = re.compile(
    r"sigma\s+\w+\s+\[pb\]\s+muR/muR0\s*=\s*" + FLOAT
    + r"(?:\s+muF/muF0\s*=\s*" + FLOAT + r")?\s*:\s*"
    + FLOAT + r"\s*\+\-\s*" + FLOAT
)


@dataclass
class Summary:
    value_sum: float = 0.0
    err2_sum: float = 0.0
    n_files: int = 0
    matched_files: int = 0
    skipped_files: int = 0
    matched_paths: List[str] = None
    # (muR factor, muF factor) -> [sum of values, sum of errors^2, number of files]
    variations: Dict[Tuple[float, float], List[float]] = None
    # False if no usable set of files was found for this contribution
    available: bool = True

    def __post_init__(self):
        if self.matched_paths is None:
            self.matched_paths = []
        if self.variations is None:
            self.variations = {}

    def add_variations(self, variations: List[Tuple[Tuple[float, float], float, float]]) -> None:
        for factor, value, err in variations:
            entry = self.variations.setdefault(factor, [0.0, 0.0, 0])
            entry[0] += value
            entry[1] += err * err
            entry[2] += 1

    def add(self, value: float, err: float, path: str) -> None:
        self.value_sum += value
        self.err2_sum += err * err
        self.matched_files += 1
        self.matched_paths.append(path)

    @property
    def error(self) -> float:
        return math.sqrt(self.err2_sum)


def extract_sigma(path: str) -> Tuple[float, float]:
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        text = f.read()
    match = SIGMA_RE.search(text)
    if not match:
        raise ValueError(f"No sigma line found in {path}")
    return float(match.group(1)), float(match.group(2))


def extract_sigma_variations(path: str) -> List[Tuple[Tuple[float, float], float, float]]:
    """((muR factor, muF factor), value, error) of each scale-variation line;
    muF factor 1 if the line has none (muR variations only)."""
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        text = f.read()
    return [((round(float(m.group(1)), 6),
              round(float(m.group(2)), 6) if m.group(2) is not None else 1.0),
             float(m.group(3)), float(m.group(4)))
            for m in SIGMA_VAR_RE.finditer(text)]


# ---------------------------------------------------------------------------
# Files recorded by bin/run_all.py in Events/run_manifest.json: only the
# outputs of the last complete run of each contribution are used, so that
# other files left in the P1_* directories (other seeds, saved copies, sectors
# that no longer exist) are never summed. '--all-files' uses every matching file.
EVENTS_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.dirname(EVENTS_DIR)
MANIFEST = os.path.join(EVENTS_DIR, "run_manifest.json")
TOP_PATTERNS = {"B": "LO_*", "V": "NLO_V_*", "R": "NLO_R_*"}
NAMES = {"B": "Born", "V": "virtual", "R": "real"}
# as in bin/run_all.py
MODES = {
    "B": ["B"],
    "V": ["V"],
    "R": ["R"],
    "NLOcorr": ["V", "R"],
    "ALL": ["B", "V", "R"],
}
# sums of contributions: (title in results.txt, name of the histograms, parts)
COMBINATIONS = [("B + V + R", "NLO", ("B", "V", "R")),
                ("V + R", "NLOcorr", ("V", "R"))]


def files_from_manifest(label, key):
    """Files of kind key ('results' or 'plot') of the last complete run of
    contribution label: (files, errors, warnings); files is None if they
    cannot be used."""
    try:
        with open(MANIFEST) as f:
            manifest = json.load(f)
    except (IOError, ValueError):
        return None, [f"no {os.path.relpath(MANIFEST, OUTPUT_DIR)} (written by bin/run_all.py)"], []
    entry = manifest.get("contributions", {}).get(label)
    if entry is None:
        return None, [f"no complete {NAMES[label]} run recorded by bin/run_all.py "
                      f"(the last one failed, was interrupted or never ran)"], []
    files, errors, warnings = [], [], []
    for run in entry["runs"]:
        path = os.path.join(OUTPUT_DIR, run[key])
        if not os.path.isfile(path):
            errors.append(f"missing: {run[key]}")
            continue
        if abs(os.path.getmtime(path) - run[key + "_mtime"]) > 1e-3:
            warnings.append(f"modified after the run_all.py run: {run[key]}")
        files.append(path)
    if errors:
        return None, errors, warnings
    return files, [], warnings


def p1_files(label, file_glob):
    """Files matching file_glob in the P1_* directories of label and in their
    sector_<i>_<j> sub-folders (real emission)."""
    found = []
    for sub in ("P1_*", os.path.join("P1_*", "sector_*")):
        found += glob.glob(os.path.join(OUTPUT_DIR, TOP_PATTERNS[label], "SubProcesses", sub, file_glob))
    return sorted(found)


def other_files(label, file_glob, used):
    """Files matching file_glob in the P1_* directories of label, not used."""
    found = p1_files(label, file_glob)
    used = set(os.path.abspath(f) for f in used)
    return sorted(f for f in found if os.path.abspath(f) not in used)


def select_files(label, key, file_glob, all_files):
    """Files to use for contribution label (None if unavailable), printing
    what is used and what is ignored."""
    if all_files:
        files = p1_files(label, file_glob)
        print(f"{label}: --all-files, using all {len(files)} {file_glob} files")
        return files
    files, errors, warnings = files_from_manifest(label, key)
    for w in warnings:
        print(f"{label}: WARNING: {w}")
    if files is None:
        for e in errors:
            print(f"{label}: not used: {e}")
        return None
    ignored = other_files(label, file_glob, files)
    print(f"{label}: {len(files)} {file_glob} files of the last run_all.py run"
          + (f" ({len(ignored)} other {file_glob} files ignored)" if ignored else ""))
    return files


# ---------------------------------------------------------------------------
# results.txt
def collect_sigma(label, result_logs) -> Summary:
    summary = Summary()
    if result_logs is None:
        summary.available = False
        return summary
    for path in result_logs:
        summary.n_files += 1
        try:
            value, err = extract_sigma(path)
            summary.add(value, err, path)
            summary.add_variations(extract_sigma_variations(path))
        except ValueError:
            summary.skipped_files += 1
    return summary


def format_block(name: str, summary: Summary) -> List[str]:
    lines = [f"{name}:"]
    if not summary.available:
        lines.append("  not available (see the messages above; run 'bin/run_all.py' "
                     "or use --all-files)")
        lines.append("")
        return lines
    lines.append(f"  files found      : {summary.n_files}")
    lines.append(f"  files matched    : {summary.matched_files}")
    lines.append(f"  files skipped    : {summary.skipped_files}")
    lines.append(f"  sum sigma [pb]   : {summary.value_sum:.12g}")
    lines.append(f"  combined err [pb]: {summary.error:.12g}")
    lines.extend(format_variations(summary.variations, summary.matched_files,
                                   summary.value_sum))
    lines.append("")
    return lines


def combine_variations(*summaries: Summary) -> Dict[Tuple[float, float], List[float]]:
    combined: Dict[Tuple[float, float], List[float]] = {}
    for summary in summaries:
        for factor, (value, err2, n) in summary.variations.items():
            entry = combined.setdefault(factor, [0.0, 0.0, 0])
            entry[0] += value
            entry[1] += err2
            entry[2] += n
    return combined


def format_variations(variations: Dict[Tuple[float, float], List[float]], n_expected: int,
                      vegas_central: float) -> List[str]:
    """Values at each (muR, muF) pair of factors and, as in MG5_aMC
    (print_summary), the envelope of all of them relative to the central one:
    the reweighted value at muR/muR0 = muF/muF0 = 1 (MG5's scale[0]), or the
    vegas central value if (1, 1) is not among the factors. Without muF
    variations (all muF factors 1) only muR is shown, as before."""
    if not variations:
        return []
    lines = []
    with_muf = any(f != 1.0 for (_, f) in variations)
    # same order as the weights of the HwU files: muF outer, muR inner
    for factor in sorted(variations, key=lambda rf: (rf[1], rf[0])):
        value, err2, n = variations[factor]
        if with_muf:
            label = f"muR/muR0 = {factor[0]:g}, muF/muF0 = {factor[1]:g}"
        else:
            label = f"muR/muR0 = {factor[0]:g}"
        line = f"  {label}: sigma [pb] = {value:.12g} +- {math.sqrt(err2):.12g}"
        if n != n_expected:
            line += f"   (WARNING: only {n} of {n_expected} files have this scale factor)"
        lines.append(line)
    values = [v[0] for v in variations.values()]
    central = variations[(1.0, 1.0)][0] if (1.0, 1.0) in variations else vegas_central
    if central != 0.0:
        s_max = (max(values) / central - 1.0) * 100.0
        s_min = (1.0 - min(values) / central) * 100.0
        lines.append(f"  Scale variation (envelope of {len(values)} values): "
                     f"{central:8.3e} pb  +{s_max:0.1f}% -{s_min:0.1f}%")
    return lines


def results_lines(labels, summaries, show_files) -> List[str]:
    lines = []
    for label in labels:
        lines.extend(format_block(label, summaries[label]))
    for title, _, parts in COMBINATIONS:
        if not all(p in labels for p in parts):
            continue
        lines.append(f"{title}:")
        if all(summaries[p].available for p in parts):
            value = sum(summaries[p].value_sum for p in parts)
            err = math.sqrt(sum(summaries[p].err2_sum for p in parts))
            lines.append(f"  total sigma [pb] : {value:.12g}")
            lines.append(f"  combined err [pb]: {err:.12g}")
            lines.extend(format_variations(
                combine_variations(*(summaries[p] for p in parts)),
                sum(summaries[p].matched_files for p in parts), value))
        else:
            lines.append("  not available (missing contributions)")
        lines.append("")
    if show_files:
        for label in labels:
            lines.append(f"Matched files for {label}:")
            for path in summaries[label].matched_paths:
                lines.append(f"  {os.path.relpath(path, OUTPUT_DIR)}")
            lines.append("")
    return lines


# ---------------------------------------------------------------------------
# tests.txt: the parsing and the checks are those of bin/run_all.py
# --test_only (imported from there)
def tests_lines(labels, tol_poles, tol_limits):
    """(lines, error): summary of the test logs of the V and R contributions
    among labels."""
    import importlib.util
    path = os.path.join(OUTPUT_DIR, "bin", "run_all.py")
    try:
        spec = importlib.util.spec_from_file_location("madnklo_run_all", path)
        run_all = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(run_all)
        return run_all.collect_test_logs(tol_poles, tol_limits, labels=labels), None
    except Exception as e:
        return None, f"cannot use {os.path.relpath(path, OUTPUT_DIR)} ({e})"


def without_details(lines):
    """The summary part of the test lines (the approach to the limits, table
    by table, is only written to tests.txt)."""
    for i, line in enumerate(lines):
        if line.startswith("Details: approach to the limits"):
            return lines[:max(0, i - 2)]
    return lines


# ---------------------------------------------------------------------------
# histograms (madgraph/various/histograms.py)
HISTOGRAMS = os.path.join(OUTPUT_DIR, "..", "madgraph", "various", "histograms.py")
PLOT_GLOBS = {"B": "seed_*_plot_B.dat", "V": "seed_*_plot_V.dat", "R": "seed_*_plot_R_*_*.dat"}


def run_command(cmd, log):
    """Run cmd (in the current directory), its output appended to log."""
    with open(log, "a") as f:
        f.write("$ " + " ".join(cmd) + "\n")
        f.flush()
        subprocess.run(cmd, check=True, stdout=f, stderr=subprocess.STDOUT)


def run_histograms_sum(input_files, out_prefix, log):
    run_command([
        HISTOGRAMS,
        *input_files,
        "--sum",
        "--HwU",
        # keep the scale-variation weights ('dyn= .. muR= .. muF= ..'), which
        # histograms.py drops by default; MG5_aMC's combine_plots_HwU does the same
        "--keep_all_weights",
        "--no_open",
        f"--out={out_prefix}",
    ], log)


def run_histograms_plot(prefix, log):
    """Gnuplot output of <prefix>.HwU with the scale variations: the extra
    weights are kept (as in MG5_aMC's combine_plots_HwU) and the scale
    variation is drawn as a band in the main panel and in the lower panel of
    relative uncertainties (histograms.py default; MG5_aMC's '--band=[]' draws
    it only as lines in the lower panel).
    Writes <prefix>.gnuplot (and rewrites <prefix>.HwU with the same content),
    then runs gnuplot, if available, to produce the .ps/.pdf files.
    Returns a note on the gnuplot step."""
    run_command([
        HISTOGRAMS,
        f"{prefix}.HwU",
        "--gnuplot",
        "--keep_all_weights",
        "--no_open",
        f"--out={prefix}",
    ], log)
    if not shutil.which("gnuplot"):
        return f"gnuplot not found: run 'gnuplot {prefix}.gnuplot' to produce the plots"
    try:
        run_command(["gnuplot", f"{prefix}.gnuplot"], log)
    except subprocess.CalledProcessError:
        return f"gnuplot failed on {prefix}.gnuplot (see histograms.log)"
    return None


def make_plots(labels, plot_files):
    """Summed histograms (in the current directory) and plots of the top
    combination(s) of the mode. Returns (created, problems)."""
    log = "histograms.log"
    created, problems = [], []
    if not os.path.isfile(HISTOGRAMS):
        return created, [f"histograms.py not found at {os.path.normpath(HISTOGRAMS)}"]
    available = {}
    try:
        for label in labels:
            if plot_files[label]:
                run_histograms_sum(plot_files[label], label, log)
                created.append(f"{label}.HwU")
                available[label] = True
            else:
                problems.append(f"{label}.HwU not produced: no {NAMES[label]} histograms")
        to_plot = []
        for _, name, parts in reversed(COMBINATIONS):    # V + R, then B + V + R
            if not all(p in labels for p in parts):
                continue
            if all(available.get(p) for p in parts):
                run_histograms_sum([os.path.abspath(f"{p}.HwU") for p in parts], name, log)
                created.append(f"{name}.HwU")
                to_plot.append(name)
            else:
                problems.append(f"{name}.HwU not produced: missing contributions")
        if len(labels) == 1 and available.get(labels[0]):
            to_plot.append(labels[0])
        for name in to_plot:
            note = run_histograms_plot(name, log)
            created += [f for f in (f"{name}.gnuplot", f"{name}.pdf", f"{name}.ps")
                        if os.path.isfile(f)]
            if note:
                problems.append(note)
    except subprocess.CalledProcessError as e:
        out = [a for a in e.cmd if str(a).startswith("--out=")]
        problems.append(f"{os.path.basename(e.cmd[0])} failed (exit status {e.returncode})"
                        + (f" for {out[0][6:]}" if out else "") + "; see histograms.log")
    except OSError as e:
        problems.append(f"cannot run histograms.py ({e})")
    return created, problems


# ---------------------------------------------------------------------------
# the output directory Events/run_<MODE>_<N>/
def new_run_dir(stem):
    """Create and return Events/<stem>_<N>, N = 1, 2, ... the first free number
    (mkdir reserves the name: two collections never share a directory)."""
    i = 1
    while True:
        path = os.path.join(EVENTS_DIR, f"{stem}_{i}")
        try:
            os.mkdir(path)
            return path
        except FileExistsError:
            i += 1


def named_run_dir(name, force):
    """Events/<name> (--out); an existing one is replaced only with --force."""
    if not name or os.path.basename(name) != name or name in (".", ".."):
        sys.exit(f"Error: --out takes a directory name (created in Events/), not '{name}'")
    path = os.path.join(EVENTS_DIR, name)
    if os.path.exists(path):
        if not force:
            sys.exit(f"Error: 'Events/{name}' already exists: choose another --out name, "
                     f"or add --force to replace it.")
        if os.path.isdir(path) and not os.path.islink(path):
            shutil.rmtree(path)
        else:
            os.remove(path)
    os.mkdir(path)
    return path


def provenance_lines(labels, all_files):
    """Which run_all.py runs the collected files come from."""
    lines = [f"collected on {time.strftime('%Y-%m-%d %H:%M:%S')}"]
    if all_files:
        lines.append("--all-files: every matching file in the P1_* directories was used")
        return lines
    try:
        with open(MANIFEST) as f:
            contributions = json.load(f).get("contributions", {})
    except (IOError, ValueError):
        contributions = {}
    for label in labels:
        entry = contributions.get(label)
        when = entry["finished"] if entry else "no complete run"
        lines.append(f"{label}: bin/run_all.py run finished {when}")
    return lines


ANALYSIS_RE = re.compile(r"^\s*(\S+)\s*=\s*analysis\s*(?:[!#].*)?$", re.I)


def card_analysis():
    """The analysis file named in Cards/run_card.dat (default analysis.f), as
    in bin/run_all.py."""
    try:
        with open(os.path.join(OUTPUT_DIR, "Cards", "run_card.dat"), errors="replace") as f:
            for line in f:
                if line.lstrip().startswith("#"):
                    continue
                m = ANALYSIS_RE.match(line.rstrip("\n"))
                if m:
                    return m.group(1).strip("'\"")
    except IOError:
        pass
    return "analysis.f"


def shown(path):
    return os.path.relpath(os.path.abspath(path), OUTPUT_DIR)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Collect the cross sections, test logs and histograms of the "
                    "last bin/run_all.py run of each contribution of MODE into "
                    "Events/run_<MODE>_<N>/")
    parser.add_argument(
        "mode", nargs="?", default="ALL", choices=list(MODES),
        help="contributions to collect: B, V, R, NLOcorr = V+R, ALL = B+V+R (default ALL)")
    parser.add_argument(
        "--out", default=None,
        help="Name of the output directory (in Events/) instead of run_<MODE>_<N>")
    parser.add_argument(
        "--force", action="store_true",
        help="Replace the --out directory if it already exists")
    parser.add_argument(
        "--all-files", action="store_true",
        help="Use every results_*.log / seed_*_plot_*.dat file found, instead of only "
             "those of the last complete bin/run_all.py run")
    parser.add_argument(
        "--show-files", action="store_true",
        help="Also list the results_*.log files used in results.txt")
    parser.add_argument(
        "--no-tests", action="store_true",
        help="Do not write tests.txt (summary of the V+I pole and IR-limit test logs)")
    parser.add_argument(
        "--no-plots", action="store_true",
        help="Do not sum the histograms nor make the plots")
    parser.add_argument(
        "--tests-only", "--tests", dest="tests_only", action="store_true",
        help="Only write tests.txt")
    parser.add_argument(
        "--tol-poles", type=float, default=1e-5,
        help="tests: tolerance on |V+I|/max(|V|,|I|) for the poles (default 1e-5)")
    parser.add_argument(
        "--tol-limits", type=float, default=1e-3,
        help="tests: tolerance on the best |R-K|/|K| reached in each limit (default 1e-3)")
    args = parser.parse_args()
    labels = MODES[args.mode]
    test_labels = [label for label in labels if label in ("V", "R")]
    do_tests = bool(test_labels) and not args.no_tests
    do_results = not args.tests_only
    do_plots = do_results and not args.no_plots
    if args.tests_only and not test_labels:
        sys.exit(f"Nothing to test in mode {args.mode} (tests exist for V and R only)")

    # files first (so that an unusable set stops before creating anything)
    result_logs, plot_files = {}, {}
    if do_results:
        for label in labels:
            result_logs[label] = select_files(label, "results", "results_*.log", args.all_files)
        if do_plots:
            for label in labels:
                plot_files[label] = select_files(label, "plot", PLOT_GLOBS[label], args.all_files)
        if not any(result_logs.values()):
            sys.exit(f"\nNothing to collect for mode {args.mode}: run 'bin/run_all.py {args.mode}' "
                     f"first (or use --all-files).")

    out_dir = (named_run_dir(args.out, args.force) if args.out
               else new_run_dir(f"run_{args.mode}" + ("_tests" if args.tests_only else "")))
    os.chdir(out_dir)
    header = [f"bin/run_all.py {args.mode}: Events/{os.path.basename(out_dir)}"]
    header += provenance_lines(labels, args.all_files)
    analysis = card_analysis()
    header += [f"analysis: Common_Files/{analysis}", ""]
    created, problems = [], []

    # cards: the commands that generated the output (proc_card.dat; called
    # proc_card_mg5.dat in outputs generated before), the run card (compiled in
    # by run_all.py), the param card (read at run time), the damping parameters
    # (compiled in) and, for the virtuals, the MadLoop ones; and the analysis
    cards = [("proc_card.dat", "proc_card.dat"), ("proc_card_mg5.dat", "proc_card.dat"),
             ("run_card.dat", "run_card.dat"), ("param_card.dat", "param_card.dat"),
             ("damping_factors.inc", "damping_factors.inc")]
    if "V" in labels:
        cards.append(("MadLoopParams.dat", "MadLoopParams.dat"))
    sources = [(os.path.join(OUTPUT_DIR, "Cards", src), dst) for src, dst in cards]
    sources.append((os.path.join(OUTPUT_DIR, "Common_Files", analysis), analysis))
    for src, dst in sources:
        if os.path.isfile(src) and not os.path.exists(dst):
            shutil.copy(src, dst)
            created.append(dst)

    if do_results:
        summaries = {label: collect_sigma(label, result_logs[label]) for label in labels}
        lines = header + results_lines(labels, summaries, args.show_files)
        with open("results.txt", "w") as f:
            f.write("\n".join(lines) + "\n")
        created.insert(0, "results.txt")
        print("")
        print("\n".join(lines))

    if do_tests:
        lines, error = tests_lines(test_labels, args.tol_poles, args.tol_limits)
        if lines is None:
            problems.append(f"tests.txt not written: {error}")
        else:
            with open("tests.txt", "w") as f:
                f.write("\n".join(lines) + "\n")
            created.append("tests.txt")
            print("\n".join(without_details(lines)))
            print("(approach to the limits, lambda -> 0: see tests.txt)\n")

    if do_plots:
        print("Summing the histograms (histograms.py, output in histograms.log) ...")
        c, p = make_plots(labels, plot_files)
        created += c
        problems += p

    with open("sources.txt", "w") as f:
        f.write("\n".join(header))
        for label in labels:
            for kind, files in (("results", result_logs), ("histograms", plot_files)):
                if label in files:
                    f.write(f"\n{label} {kind}:\n")
                    for path in files[label] or []:
                        f.write("  " + os.path.relpath(path, OUTPUT_DIR) + "\n")
    created.append("sources.txt")
    if os.path.isfile("histograms.log"):
        created.append("histograms.log")

    for problem in problems:
        print(f"WARNING: {problem}")
    print(f"\nCollected in {shown(out_dir)}/: " + ", ".join(created))


if __name__ == "__main__":
    main()
