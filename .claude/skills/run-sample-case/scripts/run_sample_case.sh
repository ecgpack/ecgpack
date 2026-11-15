#!/bin/bash
# Runs one sample_input case of an ECGPACK energy code in a scratch copy
# under jobs/. The binary is built with build.bash if it does not exist yet.
#
# Usage:   run_sample_case.sh <CODE> <case> [nprocs]
# Example: run_sample_case.sh RG_0S basis_generation_He_1Se-01 1
#
# Environment variables (optional): MACHINE, TOOLCHAIN, CONFIG, PRECISION,
# LINALG select the binary (defaults: linux-generic, systemdefault, release,
# 8, netlib). REBUILD=1 forces a rebuild, which is needed after source edits.

set -euo pipefail

if [ $# -lt 2 ] || [ $# -gt 3 ]; then
  sed -n '2,10p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
  exit 1
fi

code=$1
case_name=$2
nprocs=${3:-1}
machine=${MACHINE:-linux-generic}
toolchain=${TOOLCHAIN:-systemdefault}
config=${CONFIG:-release}
precision=${PRECISION:-8}
linalg=${LINALG:-netlib}

root=$(cd "$(dirname "${BASH_SOURCE[0]}")" && git rev-parse --show-toplevel)
case_dir=$root/$code/sample_input/$case_name

if [ ! -f "$case_dir/inout.txt" ]; then
  echo "Error: $case_dir/inout.txt not found." >&2
  echo "This script runs energy-code cases only. For the off-diagonal codes follow the README.md files of the sample case." >&2
  exit 1
fi

# The number of particles is compiled into the binary
nparticles=$(tr -d '\r' < "$case_dir/inout.txt" | awk '$1=="PARTICLES" {print $2; exit}')
if ! [[ $nparticles =~ ^[0-9]+$ ]]; then
  echo "Error: could not read PARTICLES from $case_dir/inout.txt" >&2
  exit 1
fi

# build.bash inserts the machine name into the path for these machines only
machinedir=""
case $machine in
  puma|ocelote|elgato) machinedir=$machine/ ;;
esac
binary=$root/bin/$machinedir$toolchain/$config/${code}_N${nparticles}_P${precision}_${linalg}

jobdir=$root/jobs/${code}_${case_name}_$(date +%Y%m%d-%H%M%S)
mkdir -p "$jobdir"
cp -r "$case_dir"/. "$jobdir"/

if [ ! -x "$binary" ] || [ "${REBUILD:-0}" = 1 ]; then
  if [ "$code" = PG_1P ]; then
    echo "Error: build.bash cannot build PG_1P. Build it through PG_1P/Makefile and copy the executable to" >&2
    echo "  $binary" >&2
    exit 1
  fi
  echo "Building $(basename "$binary") (log: $jobdir/build.log) ..."
  if ! (cd "$root" && ./build.bash machine="$machine" toolchain="$toolchain" config="$config" code="$code" \
        nparticles="$nparticles" precision="$precision" linalg="$linalg" openmp=no) > "$jobdir/build.log" 2>&1 \
     || [ ! -x "$binary" ]; then
    echo "Error: the build failed. Last lines of $jobdir/build.log:" >&2
    tail -n 20 "$jobdir/build.log" >&2
    exit 1
  fi
fi

echo "Running $(basename "$binary") with $nprocs MPI process(es) in $jobdir ..."
if ! (cd "$jobdir" && mpirun -np "$nprocs" "$binary") > "$jobdir/output.log" 2>&1; then
  echo "Error: the run failed. Last lines of $jobdir/output.log:" >&2
  tail -n 20 "$jobdir/output.log" >&2
  exit 1
fi

echo
echo "Result (header of $jobdir/inout.txt):"
tr -d '\r' < "$jobdir/inout.txt" | grep -m 2 -E '^ *(BASIS_SIZE|CURRENT_ENERGY) '
echo
echo "Program output:   $jobdir/output.log"
echo "Reference values: $case_dir/README.md"
