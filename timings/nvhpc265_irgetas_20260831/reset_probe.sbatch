#!/bin/bash
#SBATCH --job-name=reset265
#SBATCH --account=hpcnc
#SBATCH --partition=H100
#SBATCH --qos=hpcnc-h100
#SBATCH --gres=gpu:h100:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=reset265_%j_%N.out

# Determine whether NVHPC 26.5's exit-709 is specifically caused by the
# explicit cudaDeviceReset after all application-owned allocations are freed.
set -uo pipefail
ROOT=${CAMPAIGN_ROOT:?set CAMPAIGN_ROOT}
OUT="$ROOT/results/reset_${SLURM_JOB_ID}"
SRC="$OUT/ecgpack"

source /etc/profile.d/modules.sh
module use /shared/opt/easybuild/modules/all
module purge
module load NVHPC/26.5-CUDA-13.2.0
mkdir -p "$OUT"
cp -r "$ROOT/ecgpack" "$SRC"
rm -rf "$SRC/RG_2P/release" "$SRC/RG_2P/debug"
sed -i -E 's/Glob_AllowedNumOfParticles=[0-9]+/Glob_AllowedNumOfParticles=7/' \
  "$SRC/RG_2P/src/wp_def_8.f90"
sed -i 's/^[[:space:]]*istat = cudaDeviceReset()/    ! reset omitted: NVHPC 26.5 teardown probe/' \
  "$SRC/RG_2P/src/gpu_backend.f90"

(cd "$SRC/RG_2P" && make release COMPILER=nvfortran MACHINE=irgetas PREC=8 \
  LINALG=netlib USE_CUDA=yes CUDA_ARCH=sm_90 CUF_MAXREG=255 EXEFILE=ecg_cuda) \
  >"$OUT/build.log" 2>&1
mkdir -p "$OUT/run"
cp "$ROOT/input/RG_2P_K600.txt" "$OUT/run/inout.txt"
(cd "$OUT/run" && env ECG_BENCH=1 ECG_GPU=1 ECG_GPU_EIG=0 \
  mpirun --bind-to none -np 1 "$SRC/RG_2P/release/ecg_cuda" >run.out 2>&1)
ec=$?
echo "EXIT=$ec"
grep -m1 'BENCH K=' "$OUT/run/run.out" || true
grep -iE 'Accelerator Fatal Error|CUDA_ERROR_CONTEXT_IS_DESTROYED|GPU FATAL' \
  "$OUT/run/run.out" || true
echo RESETPROBEDONE
