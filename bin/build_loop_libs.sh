#!/bin/bash
# Builds the loop-reduction libraries shipped with madnklo in vendor/:
#   OneLOop 3.6 (scalar integrals, needed by Ninja) and Ninja 1.1.0, and,
#   with WITH_COLLIER=1, COLLIER 1.0,
# into <madnklo>/HEPTools (lib/, include/, bin/), and points madnklo to them
# (ninja = ./HEPTools/lib, and collier = ./HEPTools/lib with WITH_COLLIER=1,
# in input/mg5_configuration.txt; the previous file is kept as
# input/mg5_configuration.txt.bak_loop_libs).
# Processes generated afterwards link Ninja besides CutTools (MadLoopParams.dat:
# MLReductionLib = 6|7|1, Ninja, COLLIER, CutTools: the ones not linked are
# skipped); existing processes have to be generated again.
#
# Floating-point traps: the madnklo makefiles compile the executables with
# -ffpe-trap=zero,overflow,underflow. Ninja and COLLIER produce harmless
# underflows in normal operation, which this trap turns into crashes:
#  - COLLIER: pp > t t~, u u~ > t t~ virtual, SIGFPE in reductionD (Linux);
#  - Ninja on macOS arm64 (where a trapped FP exception arrives as SIGILL,
#    and the gfortran backtrace then hangs on SIGBUS): pp > t t~ virtual.
# Ninja is patched here to evaluate with the traps disabled (and restored
# afterwards), so it needs no change in the madnklo flags. COLLIER (Fortran,
# called from MadLoop) is not: it is off by default, WITH_COLLIER=1 builds it,
# to be used only without 'underflow' in -ffpe-trap.
#
# usage, from the madnklo directory:
#     bash bin/build_loop_libs.sh
# optional environment variables: FC (default gfortran), CXX, CC (default
# clang++/clang on macOS, g++/gcc elsewhere: the C++ library must be the one
# madnklo links, -lc++ on macOS), JOBS (unused, kept for compatibility),
# WITH_COLLIER (1: build and use COLLIER too; default 0).
# Log of each step: HEPTools/build/<step>.log
set -e
MG=$(pwd)
if [ ! -f "$MG/vendor/ninja.tar.gz" ] || [ ! -f "$MG/input/mg5_configuration.txt" ]; then
  echo "run this script from the madnklo directory"; exit 1
fi
case "$MG" in *" "*) echo "the madnklo path contains spaces: autotools would fail"; exit 1;; esac
FC=${FC:-gfortran}
if [ "$(uname -s)" = Darwin ]; then CXX=${CXX:-clang++}; CC=${CC:-clang}; else CXX=${CXX:-g++}; CC=${CC:-gcc}; fi
JOBS=${JOBS:-4}
WITH_COLLIER=${WITH_COLLIER:-0}
P=$MG/HEPTools
B=$P/build
mkdir -p $P/lib $P/include $P/bin $B
for t in $FC $CXX $CC make ar ranlib python3; do
  command -v $t >/dev/null || { echo "missing: $t"; exit 1; }
done
echo "Fortran: $($FC --version | head -1)"
echo "C++:     $($CXX --version | head -1)"
# gfortran >= 10 rejects the argument mismatches of old Fortran 77 code
FEXTRA=""
fmaj=$($FC -dumpversion | cut -d. -f1)
if [ "${fmaj:-0}" -ge 10 ] 2>/dev/null; then FEXTRA="-fallow-argument-mismatch"; fi

step() { echo "--- $1 (log: HEPTools/build/$2.log)"; }

# 1. OneLOop 3.6 -> lib/libavh_olo.a, include/oneloop/avh_olo.mod
step "OneLOop 3.6" oneloop
rm -rf $B/OneLOop-3.6
tar xzf $MG/vendor/oneloop.tar.gz -C $B
(
  cd $B/OneLOop-3.6
  python3 - "$FC" <<'PY'
import re, sys
s = open("Config").read()
s = re.sub(r"(?m)^FC *=.*$", "FC = " + sys.argv[1], s)
s = re.sub(r"(?m)^FFLAGS *=.*$", "FFLAGS = -O2 -fPIC", s)
open("Config", "w").write(s)
PY
  python3 create.py && test -f libavh_olo.a && test -f avh_olo.mod
) > $B/oneloop.log 2>&1 || { echo "OneLOop failed, see HEPTools/build/oneloop.log"; exit 1; }
mkdir -p $P/include/oneloop
cp $B/OneLOop-3.6/libavh_olo.a $P/lib/
cp $B/OneLOop-3.6/*.mod $P/include/oneloop/

# 2. Ninja 1.1.0 (double precision; with OneLOop) -> lib/libninja.a,
#    include/mninja.mod (headers in include/ninja), bin/ninja-config
step "Ninja 1.1.0" ninja
rm -rf $B/ninja-1.1.0
tar xzf $MG/vendor/ninja.tar.gz -C $B
(
  cd $B/ninja-1.1.0
  # on a case-insensitive file system (macOS) the file VERSION, found through
  # -I.., hides the C++ standard header <version> included by libc++
  mv VERSION VERSION.txt
  # Ninja runs with the floating-point traps of the calling program
  # disabled (the madnklo executables trap division by zero, overflow and
  # underflow): libc++ (macOS) divides complex numbers with logb/scalbn, and
  # logb(0) raises the division-by-zero flag, i.e. SIGILL on arm64 under a
  # trap, although the result is correct. The guard saves the environment,
  # evaluates with non-stop arithmetic and restores the environment
  # (flags included) on exit.
  python3 - <<'PY'
p = "src/ninja_wraps.cc"
s = open(p).read()
a = "#include <exception>\n"
assert s.count(a) == 1
s = s.replace(a, a + """#include <cfenv>

namespace {
  // madnklo: no floating-point traps inside Ninja (see build_loop_libs.sh)
  struct NinjaFPGuard {
    std::fenv_t env;
    NinjaFPGuard() { std::feholdexcept(&env); }
    ~NinjaFPGuard() { std::fesetenv(&env); }
  };
}
""")
a = "    TensorNumerator num(n,rank,tensor);\n"
assert s.count(a) == 1
s = s.replace(a, "    NinjaFPGuard fp_guard;\n" + a)
open(p, "w").write(s)
PY
  ./configure --prefix=$P --disable-shared --enable-static --with-pic \
    --with-avholo="-L$P/lib -lavh_olo" FCINCLUDE="-I$P/include/oneloop" \
    CXX="$CXX" CC="$CC" FC="$FC" F77="$FC" \
    CXXFLAGS="-O2" FCFLAGS="-O2 $FEXTRA" FFLAGS="-O2 $FEXTRA"
  # serial make: with -j the Fortran module mninja.mod can be needed before
  # it is built ("No rule to make target mninja.mod")
  make && make install && test -f $P/lib/libninja.a
) > $B/ninja.log 2>&1 || { echo "Ninja failed, see HEPTools/build/ninja.log"; exit 1; }
# Fortran modules of Ninja (mninja.mod, ...) directly in include/: the
# makefiles of the madnklo virtual directories look for them in
# <madnklo>/HEPTools/include (-I../../../../HEPTools/include/)
find $P/include/ninja -maxdepth 1 -name "*.mod" -exec mv {} $P/include/ \;

# 3. COLLIER 1.0 (static) -> lib/libcollier.a, include/*.mod
if [ "$WITH_COLLIER" = 1 ]; then
step "COLLIER 1.0" collier
rm -rf $B/COLLIER-1.0
tar xzf $MG/vendor/collier.tar.gz -C $B
(
  cd $B/COLLIER-1.0
  # serial: the modules must be built in the order of makefile_basic
  make -f makefile_basic FC="$FC" LINK=STATIC MODE=RELEASE \
    FFLAGS="-ffixed-line-length-132 -fno-default-integer-8 -Jmodules -Imodules -O2 -funroll-loops -fPIC $FEXTRA" \
    all && test -f libcollier.a && test -f modules/collier.mod
) > $B/collier.log 2>&1 || { echo "COLLIER failed, see HEPTools/build/collier.log"; exit 1; }
cp $B/COLLIER-1.0/libcollier.a $P/lib/
ranlib $P/lib/libcollier.a
# Fortran modules (collier.mod, ...) directly in include/, as for Ninja
cp $B/COLLIER-1.0/modules/*.mod $P/include/
fi

# 4. madnklo configuration
cfg=$MG/input/mg5_configuration.txt
[ -f $cfg.bak_loop_libs ] || cp $cfg $cfg.bak_loop_libs
python3 - "$cfg" "$WITH_COLLIER" <<'PY'
import re, sys
p = sys.argv[1]
s = open(p).read()
for tool in ("ninja", "collier"):
    pat = r"(?m)^\s*#?\s*%s\s*=.*$" % tool
    new = "%s = ./HEPTools/lib" % tool
    if tool == "collier" and sys.argv[2] != "1":
        new = "collier = None"
    if re.search(pat, s):
        s = re.sub(pat, new, s, count=1)
    else:
        s += "\n" + new + "\n"
open(p, "w").write(s)
PY
echo "--- installed in HEPTools:"
ls -l $P/lib
echo "--- input/mg5_configuration.txt:"
grep -E "^ *(ninja|collier) *=" $cfg
echo "Done. Generate the processes again to link the new libraries."
