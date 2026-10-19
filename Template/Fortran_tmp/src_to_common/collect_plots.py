#!/usr/bin/env python3

import glob
import os
import shutil
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
        # keep the scale-variation weights ('dyn= .. muR= .. muF= ..'), which
        # histograms.py drops by default; MG5_aMC's combine_plots_HwU does the same
        "--keep_all_weights",
        "--no_open",
        f"--out={out_prefix}",
    ]

    print("Running:")
    print("  " + " ".join(cmd))
    subprocess.run(cmd, check=True)


def run_histograms_plot(prefix):
    """Gnuplot output of <prefix>.HwU with the scale variations: the extra
    weights are kept (as in MG5_aMC's combine_plots_HwU) and the scale
    variation is drawn as a band in the main panel and in the lower panel of
    relative uncertainties (histograms.py default; MG5_aMC's '--band=[]' draws
    it only as lines in the lower panel).
    Writes <prefix>.gnuplot (and rewrites <prefix>.HwU with the same content),
    then runs gnuplot, if available, to produce the .ps/.pdf files."""
    cmd = [
        HISTOGRAMS,
        f"{prefix}.HwU",
        "--gnuplot",
        "--keep_all_weights",
        "--no_open",
        f"--out={prefix}",
    ]
    print("Running:")
    print("  " + " ".join(cmd))
    subprocess.run(cmd, check=True)
    if shutil.which("gnuplot"):
        result = subprocess.run(["gnuplot", f"{prefix}.gnuplot"])
        if result.returncode != 0:
            print(f"Warning: gnuplot failed on {prefix}.gnuplot", file=sys.stderr)
    else:
        print(f"gnuplot not found: run 'gnuplot {prefix}.gnuplot' to produce the plots")


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

    # Plots (gnuplot) with the scale variations
    run_histograms_plot("NLOcorr")
    run_histograms_plot("NLO")

    print("\nDone. Created:")
    print("  B.HwU")
    print("  V.HwU")
    print("  R.HwU")
    print("  NLOcorr.HwU")
    print("  NLO.HwU")
    print("  NLOcorr.gnuplot, NLO.gnuplot (and their .ps/.pdf if gnuplot is available)")


if __name__ == "__main__":
    main()
