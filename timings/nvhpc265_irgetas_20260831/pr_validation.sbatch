#!/bin/bash
#SBATCH --job-name=ecgpr265
#SBATCH --account=hpcnc
#SBATCH --partition=H100
#SBATCH --qos=hpcnc-h100
#SBATCH --gres=gpu:h100:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=ecgpr265_%j_%N.out

# Build through the public batch entry point, then run one H100/cusolver cell
# for every GPU-enabled energy code. This gates toolchain loading, the Irgetas
# architecture default, numerical execution, and clean process teardown.
set -uo pipefail
ROOT=${VALIDATION_ROOT:?set VALIDATION_ROOT}
SRC="$ROOT/ecgpack"
OUT="$ROOT/results/${SLURM_JOB_ID}"
mkdir -p "$OUT"

source /etc/profile.d/modules.sh
module use /shared/opt/easybuild/modules/all
export MODULEPATH=/shared/opt/easybuild/modules/all

cd "$SRC"
./build.bash machine=irgetas toolchain=nvhpc-26.5 config=release \
  code=RG_0S,RG_1P,RG_2D nparticles=4 precision=8 linalg=netlib cuda=yes \
  >"$OUT/build_n4.log" 2>&1
./build.bash machine=irgetas toolchain=nvhpc-26.5 config=release \
  code=RG_2P nparticles=5 precision=8 linalg=netlib cuda=yes \
  >"$OUT/build_n5.log" 2>&1

# build.bash runs as a child process, so its module environment cannot update
# this job shell. Reload the same runtime before launching the built binaries.
module purge
module load NVHPC/26.5-CUDA-13.2.0

status=0
summary="$OUT/summary.txt"
: >"$summary"
record() { echo "$*" | tee -a "$summary"; }
for spec in 'RG_0S 4' 'RG_1P 4' 'RG_2D 4' 'RG_2P 5'; do
  set -- $spec
  code=$1
  nparticles=$2
  bin="$SRC/bin/nvhpc-26.5/release/${code}_N${nparticles}_P8_netlib_cuda"
  if [[ ! -x "$bin" ]]; then
    record "BUILD code=$code status=FAIL_NO_BINARY"
    status=1
    continue
  fi
  if ! grep -q -- '-gpu=cc90' "$OUT/build_n${nparticles}.log"; then
    record "BUILD code=$code status=FAIL_NOT_CC90"
    status=1
    continue
  fi
  run="$OUT/$code"
  mkdir -p "$run"
  cp "$ROOT/input/${code}_K100.txt" "$run/inout.txt"
  (cd "$run" && env ECG_GPU=1 ECG_GPU_EIG=1 timeout 300 \
    mpirun --bind-to none -np 1 "$bin" >run.out 2>&1)
  ec=$?
  energy=$(awk '$1=="CURRENT_ENERGY" {print $2; exit}' "$run/inout.txt")
  errors=$(grep -icE 'GPU FATAL|cuSOLVER FATAL|Accelerator Fatal Error|CUDA_ERROR_|MPI_ABORT|segmentation|NaN' "$run/run.out" || true)
  cell=PASS
  [[ $ec -eq 0 && $errors -eq 0 && -n "$energy" ]] || { cell=FAIL; status=1; }
  record "RESULT code=$code status=$cell exit=$ec errors=$errors energy=${energy:-NA}"
done

record "CC90_OCCURRENCES=$(grep -h -c -- '-gpu=cc90' "$OUT"/build_n*.log | awk '{s+=$1} END{print s+0}')"
record "VALIDATION_EXIT=$status"
record ECGPR265DONE
exit "$status"
