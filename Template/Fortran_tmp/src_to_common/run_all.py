import argparse
import glob
import os
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed


def run(cmd, cwd):
    subprocess.run(cmd, shell=True, check=True, cwd=cwd)


def collect_p1_paths(pattern):
    paths = []
    for path in glob.glob(pattern):
        paths.extend(sorted(glob.glob(f"{path}/P1_*")))
    return paths


def run_parallel(paths, worker, n, pass_n=False):
    if not paths:
        print("No matching P1_* directories found.")
        return

    tasks = [(p, n) for p in paths] if pass_n else paths

    with ProcessPoolExecutor(max_workers=n) as executor:
        futures = {
            executor.submit(worker, task): (task[0] if pass_n else task)
            for task in tasks
        }
        for future in as_completed(futures):
            p = futures[future]
            try:
                future.result()
                print(f"[OK] {p}")
            except Exception as e:
                print(f"[FAIL] {p}: {e}")
                raise


# ---------------- LO ----------------
def process_lo(p):
    run("make -j2", p)
    run("make", p)
    run("./Born", p)
    return p


def run_lo(n):
    paths = collect_p1_paths("../LO_*/SubProcesses")
    run_parallel(paths, process_lo, n)


# ---------------- NLO_V ----------------
def process_nlo_v(p):
    run("make -j2", p)
    run("make", p)
    run("./virtual", p)
    return p


def run_nlo_v(n):
    paths = collect_p1_paths("../NLO_V_*/SubProcesses")
    run_parallel(paths, process_nlo_v, n)


# ---------------- NLO_R ----------------
def run_sector(task):
    p, sector = task
    run(f"./{os.path.basename(sector)}", p)
    return sector


def process_nlo_r(task):
    p, n = task
    run("make -j2", p)
    run("make", p)

    sectors = sorted(glob.glob(os.path.join(p, "sector_*")))
    if not sectors:
        return p

    with ProcessPoolExecutor(max_workers=n) as executor:
        futures = {executor.submit(run_sector, (p, s)): s for s in sectors}
        for future in as_completed(futures):
            s = futures[future]
            try:
                future.result()
                print(f"   [OK] {s}")
            except Exception as e:
                print(f"   [FAIL] {s}: {e}")
                raise

    return p


def run_nlo_r(n):
    paths = collect_p1_paths("../NLO_R_*/SubProcesses")
    run_parallel(paths, process_nlo_r, n, pass_n=True)


# ---------------- Combined modes ----------------
def run_nlocorr(n):
    print("=== Running NLO_V (Virtual) ===")
    run_nlo_v(n)
    print("=== Running NLO_R (Real) ===")
    run_nlo_r(n)


def run_all(n):
    print("=== Running LO (Born) ===")
    run_lo(n)
    print("=== Running NLO_V (Virtual) ===")
    run_nlo_v(n)
    print("=== Running NLO_R (Real) ===")
    run_nlo_r(n)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run LO/NLO computations in parallel across P1_* directories"
    )
    parser.add_argument(
        "mode",
        choices=["B", "V", "R", "NLOcorr", "ALL"],
        help="B = Born, V = Virtual, R = Real, NLOcorr = V+R, ALL = B+V+R",
    )
    parser.add_argument(
        "-n",
        type=int,
        default=os.cpu_count(),
        help="Number of parallel workers for P1_* and sector_*",
    )

    args = parser.parse_args()

    if args.mode == "B":
        run_lo(args.n)
    elif args.mode == "V":
        run_nlo_v(args.n)
    elif args.mode == "R":
        run_nlo_r(args.n)
    elif args.mode == "NLOcorr":
        run_nlocorr(args.n)
    elif args.mode == "ALL":
        run_all(args.n)
