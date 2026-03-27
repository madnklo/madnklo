#!/usr/bin/env python3

import glob
import os
import subprocess
import sys


HISTOGRAMS = os.path.join("../..", "madgraph", "various", "histograms.py")

B_PATTERN = os.path.join("../LO_*", "SubProcesses", "P1_*", "seed_*_plot_B.dat")
V_PATTERN = os.path.join("../NLO_V_*", "SubProcesses", "P1_*", "seed_*_plot_V.dat")
R_PATTERN = os.path.join("../NLO_R_*", "SubProcesses", "P1_*", "seed_*_plot_*_*.dat")


def run_histograms_sum(input_files, out_prefix):
    if not input_files:
        raise RuntimeError(f"No input files found for '{out_prefix}'")

    cmd = [
        HISTOGRAMS,
        *input_files,
        "--sum",
        "--HwU",
        "--no_open",
        f"--out={out_prefix}",
    ]

    print("Running:")
    print("  " + " ".join(cmd))
    subprocess.run(cmd, check=True)


def main():
    if not os.path.isfile(HISTOGRAMS):
        print(f"Error: could not find histograms.py at '{HISTOGRAMS}'", file=sys.stderr)
        sys.exit(1)

    b_files = sorted(glob.glob(B_PATTERN))
    v_files = sorted(glob.glob(V_PATTERN))
    r_files = sorted(glob.glob(R_PATTERN))

    print(f"Found {len(b_files)} Born files.")
    print(f"Found {len(v_files)} Virtual files.")
    print(f"Found {len(r_files)} Real files.")

    if not b_files:
        print("Error: no Born files found.", file=sys.stderr)
        sys.exit(1)
    if not v_files:
        print("Error: no Virtual files found.", file=sys.stderr)
        sys.exit(1)
    if not r_files:
        print("Error: no Real files found.", file=sys.stderr)
        sys.exit(1)

    # Separate sums
    run_histograms_sum(b_files, "B")
    run_histograms_sum(v_files, "V")
    run_histograms_sum(r_files, "R")

    # V + R
    run_histograms_sum(
        [os.path.abspath("V.HwU"), os.path.abspath("R.HwU")],
        "NLOcorr",
    )

    # B + V + R
    run_histograms_sum(
        [
            os.path.abspath("B.HwU"),
            os.path.abspath("V.HwU"),
            os.path.abspath("R.HwU"),
        ],
        "NLO",
    )

    print("\nDone. Created:")
    print("  B.HwU")
    print("  V.HwU")
    print("  R.HwU")
    print("  NLOcorr.HwU")
    print("  NLO.HwU")


if __name__ == "__main__":
    main()
