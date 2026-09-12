#!/bin/bash
#SBATCH --job-name=ecg265bench
#SBATCH --account=hpcnc
#SBATCH --partition=H100
#SBATCH --qos=hpcnc-h100
#SBATCH --gres=gpu:h100:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=ecg265bench_%j_%N.out

# NVHPC 26.5 rerun of the two substantial fixed-basis cells from the
# 2026-08-25 V100 campaign: RG_0S Carbon K=1000 and RG_2P Carbon K=600.
set -uo pipefail

ROOT=${CAMPAIGN_ROOT:?set CAMPAIGN_ROOT}
SRC="$ROOT/ecgpack"
INPUTS="$ROOT/input"
OUT="$ROOT/results/bench_${SLURM_JOB_ID}"

source /etc/profile.d/modules.sh
module use /shared/opt/easybuild/modules/all
module purge
module load NVHPC/26.5-CUDA-13.2.0
export OMP_NUM_THREADS=1
mkdir -p "$OUT"

{
  echo "source_commit=$(git -C "$SRC" rev-parse HEAD)"
  echo "node=$(hostname)"
  echo "started=$(date --iso-8601=seconds)"
  nvfortran --version | head -2
  nvcc --version | tail -1
  mpirun --version | head -2
  nvidia-smi --query-gpu=name,driver_version,memory.total,compute_cap --format=csv,noheader
} | tee "$OUT/manifest.txt"

build_one() {
  local code=$1 nparticles=$2
  local d="$SRC/$code"
  sed -i -E "s/Glob_AllowedNumOfParticles=[0-9]+/Glob_AllowedNumOfParticles=$nparticles/" \
    "$d/src/wp_def_8.f90"
  (cd "$d" && make cleanest >/dev/null 2>&1 || true)
  (cd "$d" && make release COMPILER=nvfortran MACHINE=irgetas PREC=8 \
    LINALG=netlib USE_CUDA=yes CUDA_ARCH=sm_90 CUF_MAXREG=255 \
    EXEFILE=ecg_cuda) >"$OUT/${code}_build.log" 2>&1
}

run_one() {
  local code=$1 k=$2 arm=$3 rep=$4
  local run="$OUT/${code}_${arm}_r${rep}"
  mkdir -p "$run"
  cp "$INPUTS/${code}_K${k}.txt" "$run/inout.txt"
  local t0 t1 wall ec
  t0=$(date +%s.%N)
  case "$arm" in
    cpu)     (cd "$run" && env ECG_BENCH=1 ECG_GPU=0 ECG_GPU_EIG=0 timeout 900 \
                 mpirun --bind-to none -np 1 "$SRC/$code/release/ecg_cuda" >run.out 2>&1); ec=$? ;;
    gpu_me)  (cd "$run" && env ECG_BENCH=1 ECG_GPU=1 ECG_GPU_EIG=0 timeout 900 \
                 mpirun --bind-to none -np 1 "$SRC/$code/release/ecg_cuda" >run.out 2>&1); ec=$? ;;
    gpu_all) (cd "$run" && env ECG_BENCH=1 ECG_GPU=1 ECG_GPU_EIG=1 timeout 900 \
                 mpirun --bind-to none -np 1 "$SRC/$code/release/ecg_cuda" >run.out 2>&1); ec=$? ;;
  esac
  t1=$(date +%s.%N)
  wall=$(awk -v a="$t0" -v b="$t1" 'BEGIN {printf "%.3f",b-a}')
  printf 'RESULT code=%s K=%s arm=%s rep=%s exit=%s wall=%ss ' \
    "$code" "$k" "$arm" "$rep" "$ec" "$wall"
  grep -m1 'BENCH K=' "$run/run.out" | tr -s ' ' || echo 'NO_BENCH_LINE'
}

build_one RG_0S 7
build_one RG_2P 7
for rep in 1 2 3; do
  for arm in cpu gpu_me gpu_all; do
    run_one RG_0S 1000 "$arm" "$rep"
    run_one RG_2P 600 "$arm" "$rep"
  done
done | tee "$OUT/summary.txt"

echo "finished=$(date --iso-8601=seconds)" | tee -a "$OUT/summary.txt"
echo ECG265BENCHDONE | tee -a "$OUT/summary.txt"
