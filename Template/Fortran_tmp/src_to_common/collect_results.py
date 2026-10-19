#!/usr/bin/env python3

import argparse
import glob
import json
import math
import os
import re
import sys
import time
from dataclasses import dataclass
from typing import Dict, List, Tuple


FLOAT = r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)"

SIGMA_RE = re.compile(
    r"sigma\s+\w+\s+\[pb\]\s*=\s*" + FLOAT + r"\s*\+\-\s*" + FLOAT
)

# renormalisation-scale variations (rw_rscale in run_card.dat):
#   ' sigma <tag> [pb] muR/muR0 =   2.00000 :  value +- error'
SIGMA_VAR_RE = re.compile(
    r"sigma\s+\w+\s+\[pb\]\s+muR/muR0\s*=\s*" + FLOAT + r"\s*:\s*"
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
    # scale factor -> [sum of values, sum of errors^2, number of files]
    variations: Dict[float, List[float]] = None
    # False if no usable set of files was found for this contribution
    available: bool = True

    def __post_init__(self):
        if self.matched_paths is None:
            self.matched_paths = []
        if self.variations is None:
            self.variations = {}

    def add_variations(self, variations: List[Tuple[float, float, float]]) -> None:
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

    value = float(match.group(1))
    error = float(match.group(2))
    return value, error


def extract_sigma_variations(path: str) -> List[Tuple[float, float, float]]:
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        text = f.read()
    return [(round(float(m.group(1)), 6), float(m.group(2)), float(m.group(3)))
            for m in SIGMA_VAR_RE.finditer(text)]


# ---------------------------------------------------------------------------
# Files recorded by bin/run_all.py in Events/run_manifest.json: only the
# outputs of the last complete run of each contribution are used, so that
# other files left in the P1_* directories (other seeds, saved copies, sectors
# that no longer exist) are never summed. '--all-files' restores the old
# behaviour (every matching file is used).
EVENTS_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.dirname(EVENTS_DIR)
MANIFEST = os.path.join(EVENTS_DIR, "run_manifest.json")
TOP_PATTERNS = {"B": "LO_*", "V": "NLO_V_*", "R": "NLO_R_*"}
NAMES = {"B": "Born", "V": "virtual", "R": "real"}


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


def other_files(label, file_glob, used):
    """Files matching file_glob in the P1_* directories of label, not used."""
    found = glob.glob(os.path.join(OUTPUT_DIR, TOP_PATTERNS[label], "SubProcesses", "P1_*", file_glob))
    used = set(os.path.abspath(f) for f in used)
    return sorted(f for f in found if os.path.abspath(f) not in used)


def select_files(label, key, file_glob, all_files):
    """Files to use for contribution label (None if unavailable), printing
    what is used and what is ignored."""
    if all_files:
        files = sorted(glob.glob(os.path.join(OUTPUT_DIR, TOP_PATTERNS[label],
                                              "SubProcesses", "P1_*", file_glob)))
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
    print(f"{label}: {len(files)} files of the last run_all.py run"
          + (f" ({len(ignored)} other {file_glob} files ignored)" if ignored else ""))
    return files
# ---------------------------------------------------------------------------


def collect(label: str, all_files: bool) -> Summary:
    summary = Summary()
    result_logs = select_files(label, "results", "results_*.log", all_files)
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
    lines = []
    lines.append(f"{name}:")
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


def combine_variations(*summaries: Summary) -> Dict[float, List[float]]:
    combined: Dict[float, List[float]] = {}
    for summary in summaries:
        for factor, (value, err2, n) in summary.variations.items():
            entry = combined.setdefault(factor, [0.0, 0.0, 0])
            entry[0] += value
            entry[1] += err2
            entry[2] += n
    return combined


def format_variations(variations: Dict[float, List[float]], n_expected: int,
                      vegas_central: float) -> List[str]:
    """Values at each muR factor and, as in MG5_aMC (print_summary), the
    envelope of all of them relative to the central one: the reweighted value
    at muR/muR0 = 1 (MG5's scale[0]), or the vegas central value if 1 is not
    among the factors."""
    if not variations:
        return []
    lines = []
    for factor in sorted(variations):
        value, err2, n = variations[factor]
        line = f"  muR/muR0 = {factor:g}: sigma [pb] = {value:.12g} +- {math.sqrt(err2):.12g}"
        if n != n_expected:
            line += f"   (WARNING: only {n} of {n_expected} files have this scale factor)"
        lines.append(line)
    values = [v[0] for v in variations.values()]
    central = variations[1.0][0] if 1.0 in variations else vegas_central
    if central != 0.0:
        s_max = (max(values) / central - 1.0) * 100.0
        s_min = (1.0 - min(values) / central) * 100.0
        lines.append(f"  Scale variation (envelope of {len(values)} values): "
                     f"{central:8.3e} pb  +{s_max:0.1f}% -{s_min:0.1f}%")
    return lines

def next_free_name(stem, ext=""):
    """stem_1<ext>, stem_2<ext>, ...: the first name that does not exist yet."""
    i = 1
    while os.path.exists(f"{stem}_{i}{ext}"):
        i += 1
    return f"{stem}_{i}{ext}"


def output_name(out, stem, ext, force):
    """Name of the output: --out if given (an existing one is only replaced
    with --force), otherwise the next free stem_N<ext> in Events/."""
    if out is None:
        return next_free_name(stem, ext)
    name = out if (not ext or out.endswith(ext)) else out + ext
    if os.path.exists(name) and not force:
        sys.exit(f"Error: '{name}' already exists: choose another --out name, "
                 f"or add --force to replace it.")
    return name


def provenance_lines(all_files):
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
    for label in ("B", "V", "R"):
        entry = contributions.get(label)
        when = entry["finished"] if entry else "no complete run"
        lines.append(f"{label}: bin/run_all.py run finished {when}")
    return lines


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Collect sigma values from results_*.log under "
            "LO_*, NLO_V_*, NLO_R_* and compute sums/errors."
        )
    )
    parser.add_argument(
        "--show-files",
        action="store_true",
        help="Print matched file paths for each category.",
    )
    parser.add_argument(
        "--all-files",
        action="store_true",
        help="Use every results_*.log file found, instead of only those of the "
             "last complete bin/run_all.py run (old behaviour)",
    )
    parser.add_argument(
        "--out",
        default=None,
        help="Name of the output (in Events/); default: the next free results_N.txt (N = 1, 2, ...)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Replace the --out output if it already exists",
    )
    args = parser.parse_args()
    os.chdir(EVENTS_DIR)
    out_file = output_name(args.out, "results", ".txt", args.force)

    summaries = {label: collect(label, args.all_files) for label in ("B", "V", "R")}
    print("")

    total_all_value = summaries["B"].value_sum + summaries["V"].value_sum + summaries["R"].value_sum
    total_all_err = math.sqrt(
        summaries["B"].err2_sum + summaries["V"].err2_sum + summaries["R"].err2_sum
    )

    total_nlocorr_value = summaries["V"].value_sum + summaries["R"].value_sum
    total_nlocorr_err = math.sqrt(summaries["V"].err2_sum + summaries["R"].err2_sum)

    # Collect output
    output_lines = provenance_lines(args.all_files) + [""]

    for name in ("B", "V", "R"):
        output_lines.extend(format_block(name, summaries[name]))

    output_lines.append("B + V + R:")
    if all(summaries[n].available for n in ("B", "V", "R")):
        output_lines.append(f"  total sigma [pb] : {total_all_value:.12g}")
        output_lines.append(f"  combined err [pb]: {total_all_err:.12g}")
        output_lines.extend(format_variations(
            combine_variations(summaries["B"], summaries["V"], summaries["R"]),
            sum(summaries[n].matched_files for n in ("B", "V", "R")), total_all_value))
    else:
        output_lines.append("  not available (missing contributions)")
    output_lines.append("")

    output_lines.append("V + R:")
    if all(summaries[n].available for n in ("V", "R")):
        output_lines.append(f"  total sigma [pb] : {total_nlocorr_value:.12g}")
        output_lines.append(f"  combined err [pb]: {total_nlocorr_err:.12g}")
        output_lines.extend(format_variations(
            combine_variations(summaries["V"], summaries["R"]),
            sum(summaries[n].matched_files for n in ("V", "R")), total_nlocorr_value))
    else:
        output_lines.append("  not available (missing contributions)")
    output_lines.append("")

    if args.show_files:
        for name in ("B", "V", "R"):
            output_lines.append(f"Matched files for {name}:")
            for path in summaries[name].matched_paths:
                output_lines.append(f"  {path}")
            output_lines.append("")

    # Print to screen
    for line in output_lines:
        print(line)

    # Write to file (never overwriting an older one, unless --out ... --force)
    with open(out_file, "w") as f:
        f.write("\n".join(output_lines) + "\n")

    print(f"\nResults written to Events/{out_file}")

if __name__ == "__main__":
    main()
