import argparse
import glob
import os
import subprocess
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
    run("make", p)
    run("./Born", p)
    return p


def run_lo(n):
    paths = collect_p1_paths("../LO_*/SubProcesses")
    run_parallel(paths, process_lo, n)


# ---------------- NLO_V ----------------
def process_nlo_v(p):
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

    check_onshell_widths()
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
