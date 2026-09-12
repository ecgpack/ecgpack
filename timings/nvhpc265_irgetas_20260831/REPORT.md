# NVHPC 26.5 / H100 rerun and compiler-defect gate

- Date: 2026-08-31
- Cluster/job nodes: Irgetas `g001`, jobs `23914`, `23916`, `23917`
- GPU: NVIDIA H100 80 GB HBM3, compute capability 9.0
- Driver: 610.43.02
- Toolchain: NVHPC 26.5, CUDA 13.2, HPC-X / OpenMPI 5.0.10rc2
- Benchmark source: `gpu-all-codes-bench` at `e9de0fa`
- Compiler-reproducer source: `matelem/cuda` at `5864aab`
- Build: `-O3 -tp=native -cuda -gpu=cc90,maxregcount:255,lto`

## Result

NVHPC 26.5 has **not** fixed either silent wrong-code defect previously confirmed
on 25.9 and 26.3. Keep both one-line runtime-bound source workarounds.

The H100 matrix-element path is about 3.8--3.9x faster than the prior V100 cells,
but 26.5 exposes an independent process-teardown failure when ECGPACK explicitly
calls `cudaDeviceReset()`.

## Fixed-basis benchmark rerun

Same source and input decks as the 2026-08-25 V100 campaign. Three repetitions;
medians shown. `GPU ME` uses the host eigensolver, while `GPU all` enables cuSOLVER.

| Code / basis | Phase | V100 / NVHPC 25.9 (s) | H100 / NVHPC 26.5 (s) | V100 / H100 |
|---|---|---:|---:|---:|
| RG_0S Carbon K=1000 | GPU ME matrix | 0.9235 | 0.2403 | **3.84x** |
| RG_0S Carbon K=1000 | GPU-all matrix | 0.9157 | 0.2400 | **3.82x** |
| RG_0S Carbon K=1000 | cuSOLVER | 0.1004 | 0.1332 | **0.75x** |
| RG_2P Carbon K=600 | GPU ME matrix | 0.9107 | 0.2329 | **3.91x** |
| RG_2P Carbon K=600 | GPU-all matrix | 0.9091 | 0.2318 | **3.92x** |
| RG_2P Carbon K=600 | cuSOLVER | 0.0948 | 0.1294 | **0.73x** |

The comparison changes GPU architecture, compiler, CUDA and driver together; it
does not isolate the compiler's contribution. The H100 matrix result is robust,
but the cuSOLVER regression should be profiled separately before enabling
`ECG_GPU_EIG=1` by default.

Input SHA-256:

```text
RG_0S K=1000  22c9541f0f809424a68abdf3111bdb65a5517fb0892b2c945fd995cf318e5b7d
RG_2P K=600   5d604cd2879f015e8e323895b5be642eaae25070767acc4ae4fef0d4eb9c6d72
```

## Defect A: constant-bound host reduction, N=5

The reduced reproducer is clean at `-O1` and fails identically at `-O2` and
`-O3`, although runtime `n` and parameter `nn` are both 4:

```text
                         runtime n             constant nn
Vkl                 -5.1554360571275870    -2.1192447057525272
Hkl                 -2.5082961369337458     0.5278952144413140
maximum relative difference                    1.2105
```

The full self-check reproduces the same values recorded on 25.9/26.3:

| Code | Independent expected H | Unfixed constant-bound H | Verdict |
|---|---:|---:|---|
| RG_1P | -138885.633170 | -202693.795168 | wrong |
| RG_2D | -75811.076904 | -120974.243940 | wrong |
| RG_2P | -77557.830920 | -120972.623225 | wrong |

`MatrixElementsAll` and the H100 kernel agree with the expected column. Only
nvfortran's host `MatrixElementsHS` copy is wrong. The shipped `do j=1,n`
workaround makes all three agree.

## Defect B: RG_1P CUDA gradient, N >= 7

| Size | Source | Energy deviation | Gradient deviation | Verdict |
|---:|---|---:|---:|---|
| 7 | shipped workaround | 8.90e-15 | 8.90e-15 | pass |
| 7 | workaround reversed | 8.90e-15 | **2.48e-01** | wrong |
| 9 | shipped workaround | 2.67e-12 | 2.67e-12 | pass |
| 9 | workaround reversed | 2.67e-12 | **2.00e-01** | wrong |

The signature is unchanged: energy remains clean, only the derivative kernel is
miscompiled, and reverting one outer loop to runtime `n` clears it.

## New 26.5 teardown incompatibility

Every production ECGPACK GPU cell completed its numerical work and printed the
benchmark line, then returned exit 1:

```text
Program has stopped
Accelerator Fatal Error: call to cuMemFree returned error 709
(CUDA_ERROR_CONTEXT_IS_DESTROYED): Context is destroyed or not yet created
```

All application-owned device allocations are explicitly deallocated before
`cudaDeviceReset()`. A controlled RG_2P K=600 rebuild that omitted only the
explicit reset produced the same result and returned `EXIT=0`. Therefore the
safe compatibility change is to stop resetting the CUDA context during normal
process teardown; the OS/runtime already owns final context destruction.

## Scope not covered

This run did not repeat the older Oxygen register-cap sweep below 255. All tests
kept the correctness requirement `maxregcount:255`; do not lower it based on
these results.

## PR fix validation

Irgetas job `23922` validated the resulting PR changes from a clean checkout.
It built through `build.bash` with `machine=irgetas toolchain=nvhpc-26.5
cuda=yes`, rather than bypassing the public build entry point. All 60 CUDA
compile/link command occurrences targeted `cc90` without a `cuda_arch` override.

Each CUDA-enabled energy code then ran a K=100 cell with matrix assembly and
cuSOLVER enabled:

```text
RESULT code=RG_0S status=PASS exit=0 errors=0 energy=-0.7477376516187604E+01
RESULT code=RG_1P status=PASS exit=0 errors=0 energy=-0.7409341827026370E+01
RESULT code=RG_2D status=PASS exit=0 errors=0 energy=-0.7334818842216896E+01
RESULT code=RG_2P status=PASS exit=0 errors=0 energy=-0.1439293058919524E+02
```

No run emitted CUDA error 709 or another GPU/cuSOLVER error. This validates the
removal of the explicit reset across all four sibling backends.
