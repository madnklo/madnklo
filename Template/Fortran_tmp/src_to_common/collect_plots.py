#!/usr/bin/env python3

import argparse
import glob
import json
import os
import shutil
import subprocess
import sys
import time


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

HISTOGRAMS = os.path.join(OUTPUT_DIR, "..", "madgraph", "various", "histograms.py")
PLOT_GLOBS = {"B": "seed_*_plot_B.dat", "V": "seed_*_plot_V.dat", "R": "seed_*_plot_R_*_*.dat"}


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


def main():
    parser = argparse.ArgumentParser(
        description="Sum the histograms of the B, V and R contributions and plot them")
    parser.add_argument(
        "--all-files",
        action="store_true",
        help="Use every seed_*_plot_*.dat file found, instead of only those of the "
             "last complete bin/run_all.py run (old behaviour)",
    )
    parser.add_argument(
        "--out",
        default=None,
        help="Name of the output (in Events/); default: the next free directory plots_N (N = 1, 2, ...)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Replace the --out output if it already exists",
    )
    args = parser.parse_args()
    os.chdir(EVENTS_DIR)
    out_dir = output_name(args.out, "plots", "", args.force)

    if not os.path.isfile(HISTOGRAMS):
        print(f"Error: could not find histograms.py at '{HISTOGRAMS}'", file=sys.stderr)
        sys.exit(1)

    files = {label: select_files(label, "plot", PLOT_GLOBS[label], args.all_files)
             for label in ("B", "V", "R")}
    if not any(files.values()):
        print("Error: nothing to sum.", file=sys.stderr)
        sys.exit(1)

    # all the outputs go to a new directory (never overwriting older plots,
    # unless --out ... --force), with the list of the files summed
    os.makedirs(out_dir, exist_ok=True)
    os.chdir(out_dir)
    with open("sources.txt", "w") as f:
        f.write("\n".join(provenance_lines(args.all_files)) + "\n")
        for label in ("B", "V", "R"):
            f.write(f"\n{label}:\n")
            for path in files[label] or []:
                f.write("  " + os.path.relpath(path, OUTPUT_DIR) + "\n")
    created = []
    for label in ("B", "V", "R"):
        if files[label]:
            run_histograms_sum(files[label], label)
            created.append(f"{label}.HwU")
        else:
            print(f"{label}.HwU not produced: no {NAMES[label]} histograms to sum")

    # V + R, B + V + R
    for out, parts in (("NLOcorr", ("V", "R")), ("NLO", ("B", "V", "R"))):
        if all(files[label] for label in parts):
            run_histograms_sum([os.path.abspath(f"{label}.HwU") for label in parts], out)
            # plots (gnuplot) with the scale variations
            run_histograms_plot(out)
            created.append(f"{out}.HwU, {out}.gnuplot (and .ps/.pdf if gnuplot is available)")
        else:
            print(f"{out} not produced: missing contributions")

    print(f"\nDone. Created in Events/{out_dir}/:")
    for c in created + ["sources.txt (the files summed)"]:
        print("  " + c)


if __name__ == "__main__":
    main()
