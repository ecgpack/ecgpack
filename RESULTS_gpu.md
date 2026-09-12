# ECGPACK on a V100: production GPU-port result

The native CUDA Fortran path is implemented and cluster-validated for the four real-ECG
energy codes: `RG_0S`, `RG_1P`, `RG_2D`, and `RG_2P`. It accelerates H/S assembly with
and without nonlinear-parameter derivatives and optionally routes generalized symmetric
eigensolves through cuSOLVER. CUDA-enabled binaries remain CPU-only unless selected with
`ECG_GPU=1` or `ECG_GPU_EIG=1`.

Source tested: `gpu-all-codes-bench` at `e9de0fa`, which is the production implementation
at `b8cf801` plus timing/full-precision-output instrumentation. Detailed provenance, raw
repetitions, and correctness deltas are in
[`timings/gpu_port_gn04_8_nvhpc-25.9.md`](timings/gpu_port_gn04_8_nvhpc-25.9.md).

## One V100 against one 64-core CPU node

Machine: shabyt `gn04`, 2 x AMD EPYC 7452, Tesla V100-PCIE-32GB. CPU binaries use
gfortran 14.3 + OpenBLAS; GPU binaries use NVHPC 25.9 + CUDA 12.9.1. Five trials,
median reported. The physical stored bases are deliberately representative rather than
equal-work synthetic inputs.

| Code / case | CPU64 H/S | V100 H/S | matrix speedup | CPU64 eig | cuSOLVER eig | wall CPU64 / GPU |
|---|---:|---:|---:|---:|---:|---:|
| `RG_0S` Carbon, K=1000 | 4.4978 s | 0.9235 s | **4.87x** | 0.1242 s | 0.1004 s | 9.336 / 3.510 s |
| `RG_1P` Lithium, K=250 | 0.0074 s | 0.0045 s | **1.64x** | 0.0075 s | 0.0715 s | 4.477 / 2.138 s* |
| `RG_2D` Lithium, K=100 | 0.0014 s | 0.0039 s | **0.36x** | 0.0034 s | 0.0711 s | 4.578 / 2.624 s* |
| `RG_2P` Carbon, K=600 | 2.4797 s | 0.9107 s | **2.72x** | 0.0364 s | 0.0948 s | 7.053 / 3.459 s |

`*` For the small cases, process startup and I/O dominate wall time; the H/S phase is the
meaningful comparison.

## What this says

A V100 is materially faster than the full 64-core node for the matrix-dominated Carbon
cases: 4.87x for `RG_0S` K=1000 and 2.72x for `RG_2P` K=600. `RG_1P` K=250 is just
above crossover. `RG_2D` K=100 is below it and is 2.8x faster on the CPU node.

The relevant scaling variable is not K alone. Full H/S work grows with K(K+1)/2 and with
the number of symmetry-projector terms. The Carbon cases expose enough independent work to
fill the V100; the small Lithium `RG_2D` sample does not.

cuSOLVER reaches parity only at K around 1000 in this dataset. It is 9.5x slower than the
CPU solve at K=250, 20.9x slower at K=100, and 2.6x slower at K=600. Keeping
`ECG_GPU_EIG=1` opt-in is the correct default.

## Derivative path

The repeated method-I `OPT_CYCLE` gate used identical bounded work in every cell: 89
energy builds and 8 gradient builds. Three trials, median milliseconds per build.

| Code / case | CPU8 energy / GPU | speedup | CPU8 gradient / GPU | speedup |
|---|---:|---:|---:|---:|
| `RG_0S` K=1000 | 235.180 / 13.079 | **17.98x** | 77.750 / 10.875 | **7.15x** |
| `RG_1P` K=250 | 0.112 / 0.157 | 0.71x | 0.125 / 0.125 | 1.00x* |
| `RG_2D` K=100 | 0.045 / 0.135 | 0.33x | 0.125 / 0.250 | 0.50x* |
| `RG_2P` K=600 | 185.472 / 13.382 | **13.86x** | 96.500 / 14.750 | **6.54x** |

`*` These sub-millisecond clocks are at timer resolution. They establish only that the
small bases are below crossover. A current-source CPU64 derivative sweep was not run.

## Correctness and reproducibility

All **144 valid runtime cells passed**. Maximum absolute energy differences for CPU64
versus GPU matrix assembly were:

| `RG_0S` | `RG_1P` | `RG_2D` | `RG_2P` |
|---:|---:|---:|---:|
| 4.460e-10 | 1.776e-15 | 0 | 8.527e-14 |

With cuSOLVER included, the maxima were 7.585e-10, 4.232e-12, 7.629e-13, and
1.631e-11 respectively. These are consistent with expected floating-point reassociation.

Default GPU matrix and gradient kernels accumulate symmetry terms with `atomicAdd`.
Repeated Carbon runs therefore agree scientifically but are not bit-reproducible.
`ECG_DETERM=1` made the fixed-basis H/S energy string identical in 3/3 runs for every
code, but it does not cover gradients and requires `16 * NumYHYTerms` bytes of dynamic
shared memory per block.

Optimization endpoints are self-consistency evidence rather than a direct gradient proof.
An analytic-versus-finite-difference CUDA gradient test remains open.

## Scope of the production claim

Supported and gated:

- double precision (`PREC=8`), NVHPC/CUDA;
- `RG_0S`, `RG_1P`, `RG_2D`, `RG_2P`;
- fixed and derivative H/S assembly;
- CPU fallback in the same CUDA-enabled binary;
- optional cuSOLVER eigensolves;
- one rank per visible GPU.

Not yet established:

- larger physical `RG_1P` and especially `RG_2D` crossover points;
- multi-rank/multi-GPU performance and placement;
- direct finite-difference gradient validation;
- GPU ports for `CG_0S` or off-diagonal matrix-element codes.

The upstreamable branch is `codex/gpu-all-codes-ready`: the eight production commits,
documentation corrections, and these validation reports, without benchmark instrumentation.
