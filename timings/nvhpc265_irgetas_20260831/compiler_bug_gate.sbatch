#!/bin/bash
#SBATCH --job-name=nvbug265
#SBATCH --account=hpcnc
#SBATCH --partition=H100
#SBATCH --qos=hpcnc-h100
#SBATCH --gres=gpu:h100:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=nvbug265_%j_%N.out

# Compare the shipped source workarounds with the exact source obtained by
# reversing commit 3b08f4c's two nvfortran fixes. The executable independently
# compares HS, All, and CUDA implementations at every basis/particle count.
set -uo pipefail

ROOT=${CAMPAIGN_ROOT:?set CAMPAIGN_ROOT}
BASE="$ROOT/matelem"
OUT="$ROOT/results/bugs_${SLURM_JOB_ID}"

source /etc/profile.d/modules.sh
module use /shared/opt/easybuild/modules/all
module purge
module load NVHPC/26.5-CUDA-13.2.0
export OMP_NUM_THREADS=1
mkdir -p "$OUT"

{
  echo "source_commit=$(git -C "$BASE" rev-parse HEAD)"
  echo "node=$(hostname)"
  echo "started=$(date --iso-8601=seconds)"
  nvfortran --version | head -2
  nvidia-smi --query-gpu=name,driver_version,memory.total,compute_cap --format=csv,noheader
} | tee "$OUT/manifest.txt"

# Defect A has a reduced host-only reproducer. O1 is the clean control;
# O2/O3 were wrong on 25.9 and 26.3.
(cd "$BASE/nvbug-repro-reduced" && \
  OPTS='-O1 -O2 -O3' ./build_and_run.sh) >"$OUT/reduced_A.out" 2>&1
echo "REDUCED_A_EXIT=$?" | tee "$OUT/reduced_A.status"

rm -rf "$OUT/fixed" "$OUT/unfixed"
cp -r "$BASE/branch_cuda" "$OUT/fixed"
cp -r "$BASE/branch_cuda" "$OUT/unfixed"
rm -rf "$OUT/fixed/release" "$OUT/fixed/debug" "$OUT/unfixed/release" "$OUT/unfixed/debug"
sed -i 's/TargetTime=12.0d0/TargetTime=0.05d0/' "$OUT/fixed/src/main.f90" "$OUT/unfixed/src/main.f90"
(cd "$OUT/unfixed" && git apply -R "$ROOT/input/nvfortran_two_fixes.patch")

run_variant() {
  local variant=$1
  local nparticles=$2
  local d="$OUT/$variant"
  rm -rf "$d/release"
  if ! (cd "$d" && make release COMPILER=nvfortran CUDA=yes CUDA_ARCH=sm_90 \
      PREC=8 NPARTICLES="$nparticles" EXEFILE=x) >"$OUT/build_${variant}_N${nparticles}.log" 2>&1; then
    echo "BUILD variant=$variant N=$nparticles status=FAIL"
    tail -20 "$OUT/build_${variant}_N${nparticles}.log"
    return
  fi
  echo "=== variant=$variant N=$nparticles ==="
  (cd "$d" && timeout 300 mpirun --bind-to none -np 1 ./release/x) 2>&1 \
    | grep -E '^ VALUES RG_|^ RG_[A-Za-z0-9_]+_gpu '
}

for n in 5 7 9; do
  run_variant fixed "$n"
  run_variant unfixed "$n"
done | tee "$OUT/full_gate.out"

echo "finished=$(date --iso-8601=seconds)" | tee "$OUT/summary.txt"
echo NVBUG265DONE | tee -a "$OUT/summary.txt"
