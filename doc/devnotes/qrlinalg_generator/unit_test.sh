#!/bin/bash
# Build module qrlinalg and its standalone MPI unit test (test_qrlinalg.f90)
# and run the unit suite at several process counts.
#
#   ./unit_test.sh [PREC ...]          PREC: 8, 10, 16; default 8
#
# Settings (environment variables):
#   MODE        release (default) or debug (bounds checks, trapping of invalid,
#               zero-divide, and overflow exceptions, NaN initialization)
#   NPROCS      process counts, default "1 2 3". The first one writes the
#               reference trace; every other count is compared with it. With
#               gfortran, precision 16 runs with one process only.
#   TEST_ARGS   further arguments of the driver, for example "quick",
#               "factor=serial", "factor=distributed", "real", "complex"
#   CODE_SRC    source directory that supplies qrlinalg.f90, wp_def_<PREC>.f90,
#               X1MACH.f90, BLAS.f, and LAPACK.f; default RG_0S/src
#   BUILD_ROOT  where to build and run; default jobs/qrlinalg_unit in the
#               repository root (jobs/ is not tracked by git)
#   FC          MPI Fortran compiler wrapper, default mpif90 (gfortran flags
#               are used)
#   MPIRUN      launcher, default "mpirun"; for example "mpirun --bind-to none"
#               on a node shared with other jobs
#
# Products: $BUILD_ROOT/p<PREC>[_dbg]/test_qrlinalg, unit_np<N>.log, and
# trace_np<N>.txt. The exit status is zero when every run printed
# "RESULT: ALL TESTS PASSED".
#
# Release builds use -march=native: build and run on the same kind of machine.
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
CODE_SRC="${CODE_SRC:-$REPO/RG_0S/src}"
BUILD_ROOT="${BUILD_ROOT:-$REPO/jobs/qrlinalg_unit}"
MODE="${MODE:-release}"
NPROCS="${NPROCS:-1 2 3}"
TEST_ARGS="${TEST_ARGS:-}"
FC="${FC:-mpif90}"
MPIRUN="${MPIRUN:-mpirun}"
PRECS="$*"
[ -z "$PRECS" ] && PRECS="8"
export GFORTRAN_UNBUFFERED_ALL=1

failed=0
for P in $PRECS; do
  if [ "$MODE" = debug ]; then
    B="$BUILD_ROOT/p${P}_dbg"
    OPT="-O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow -finit-real=nan -fimplicit-none"
  else
    B="$BUILD_ROOT/p${P}"
    OPT="-O3 -march=native"
  fi
  mkdir -p "$B" || exit 1
  # QRLINALG_TESTING makes the components of the state types public for the
  # white-box factorization tests of the driver. It must stay undefined in
  # an ECGPACK build. The sources are compiled from where they are.
  FCFLAGS="-J$B -I$B -cpp -DQRLINALG_TESTING $OPT"
  FIXED="-ffixed-line-length-none"
  (
    set -e
    $FC $FCFLAGS -c "$CODE_SRC/wp_def_$P.f90" -o "$B/wp_def_$P.o"
    $FC $FCFLAGS -c "$CODE_SRC/X1MACH.f90" -o "$B/X1MACH.o"
    $FC $FCFLAGS $FIXED -DQRLINALG_WP=$P -c "$CODE_SRC/BLAS.f" -o "$B/BLAS.o"
    $FC $FCFLAGS $FIXED -c "$CODE_SRC/LAPACK.f" -o "$B/LAPACK.o"
    $FC $FCFLAGS -Wall -Wno-unused-dummy-argument -Wno-maybe-uninitialized \
        -c "$CODE_SRC/qrlinalg.f90" -o "$B/qrlinalg.o"
    $FC $FCFLAGS -Wall -Wno-maybe-uninitialized \
        -c "$HERE/test_qrlinalg.f90" -o "$B/test_qrlinalg.o"
    $FC $OPT -o "$B/test_qrlinalg" "$B/test_qrlinalg.o" "$B/qrlinalg.o" \
        "$B/LAPACK.o" "$B/BLAS.o" "$B/X1MACH.o" "$B/wp_def_$P.o"
  ) > "$B/build.log" 2>&1 || { echo "precision $P ($MODE): BUILD FAILED, see $B/build.log"; failed=1; continue; }
  echo "precision $P ($MODE): built $B/test_qrlinalg"

  ref=""
  for NP in $NPROCS; do
    trace="trace_np$NP.txt"
    log="unit_np$NP.log"
    start=$(date +%s)
    if [ -z "$ref" ]; then
      (cd "$B" && $MPIRUN -np "$NP" ./test_qrlinalg unit trace=$trace $TEST_ARGS > "$log" 2>&1)
      ref=$trace
    else
      (cd "$B" && $MPIRUN -np "$NP" ./test_qrlinalg unit trace=$trace ref=$ref $TEST_ARGS > "$log" 2>&1)
    fi
    rc=$?
    result=$(grep -a 'RESULT' "$B/$log" | tail -1)
    echo "precision $P ($MODE) np=$NP exit=$rc seconds=$(( $(date +%s) - start )) : ${result:-no RESULT line, see $B/$log}"
    case "$result" in *"ALL TESTS PASSED"*) [ $rc -eq 0 ] || failed=1;; *) failed=1;; esac
  done
done
exit $failed
