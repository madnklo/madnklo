#!/usr/bin/env python3

import argparse
import glob
import math
import os
import re
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


def collect_from_pattern(top_pattern: str) -> Summary:
    summary = Summary()

    top_dirs = sorted(glob.glob(top_pattern))
    for top_dir in top_dirs:
        result_logs = sorted(glob.glob(os.path.join(top_dir, "SubProcesses", "P1_*", "results_*.log")))
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
    args = parser.parse_args()

    summaries = {
        "B": collect_from_pattern("../LO_*"),
        "V": collect_from_pattern("../NLO_V_*"),
        "R": collect_from_pattern("../NLO_R_*"),
    }

    total_all_value = summaries["B"].value_sum + summaries["V"].value_sum + summaries["R"].value_sum
    total_all_err = math.sqrt(
        summaries["B"].err2_sum + summaries["V"].err2_sum + summaries["R"].err2_sum
    )

    total_nlocorr_value = summaries["V"].value_sum + summaries["R"].value_sum
    total_nlocorr_err = math.sqrt(summaries["V"].err2_sum + summaries["R"].err2_sum)

    # Collect output
    output_lines = []

    for name in ("B", "V", "R"):
        output_lines.extend(format_block(name, summaries[name]))

    output_lines.append("B + V + R:")
    output_lines.append(f"  total sigma [pb] : {total_all_value:.12g}")
    output_lines.append(f"  combined err [pb]: {total_all_err:.12g}")
    output_lines.extend(format_variations(
        combine_variations(summaries["B"], summaries["V"], summaries["R"]),
        sum(summaries[n].matched_files for n in ("B", "V", "R")), total_all_value))
    output_lines.append("")

    output_lines.append("V + R:")
    output_lines.append(f"  total sigma [pb] : {total_nlocorr_value:.12g}")
    output_lines.append(f"  combined err [pb]: {total_nlocorr_err:.12g}")
    output_lines.extend(format_variations(
        combine_variations(summaries["V"], summaries["R"]),
        sum(summaries[n].matched_files for n in ("V", "R")), total_nlocorr_value))
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

    # Write to file
    with open("results.txt", "w") as f:
        f.write("\n".join(output_lines))

    print("\nResults written to results.txt")

if __name__ == "__main__":
    main()
