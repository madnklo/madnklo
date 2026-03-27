#!/usr/bin/env python3

import argparse
import glob
import math
import os
import re
from dataclasses import dataclass
from typing import List, Tuple


SIGMA_RE = re.compile(
    r"sigma\s+\w+\s+\[pb\]\s*=\s*"
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)\s*"
    r"\+\-\s*"
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)"
)


@dataclass
class Summary:
    value_sum: float = 0.0
    err2_sum: float = 0.0
    n_files: int = 0
    matched_files: int = 0
    skipped_files: int = 0
    matched_paths: List[str] = None

    def __post_init__(self):
        if self.matched_paths is None:
            self.matched_paths = []

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


def collect_from_pattern(top_pattern: str) -> Summary:
    summary = Summary()

    top_dirs = sorted(glob.glob(top_pattern))
    for top_dir in top_dirs:
        result_logs = sorted(glob.glob(os.path.join(top_dir, "SubProcesses", "P1_*", "results*.log")))
        for path in result_logs:
            summary.n_files += 1
            try:
                value, err = extract_sigma(path)
                summary.add(value, err, path)
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
    lines.append("")
    return lines

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Collect sigma values from results*.log under "
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
    output_lines.append("")

    output_lines.append("V + R:")
    output_lines.append(f"  total sigma [pb] : {total_nlocorr_value:.12g}")
    output_lines.append(f"  combined err [pb]: {total_nlocorr_err:.12g}")
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
