# ECGPACK production GPU-port timings

- date: 2026-08-25
- host: `gn04`, shabyt
- CPU: 2 x AMD EPYC 7452, 64 physical cores
- GPU: Tesla V100-PCIE-32GB
- precision: 8 (IEEE double precision)
- source under test: `gpu-all-codes-bench` at `e9de0fa`
- production source represented: `gpu-all-codes` at `b8cf801`
- GPU toolchain: NVHPC 25.9, CUDA 12.9.1, driver 580.178.04
- CPU toolchain: GCC/gfortran 14.3.0, OpenMPI 5.0.8, OpenBLAS from `foss/2025b`
- fixed-basis trials: 5, median reported
- derivative-path trials: 3, median reported

The benchmark instrumentation in `e9de0fa` adds clocks and full-precision output only;
it does not change the production matrix-element or eigensolver algorithms in `b8cf801`.

## What the clocks measure

- `MEE`: complete energy-path H/S assembly, including the matrix-element kernel,
  host/device transfers, MPI reduction, host orchestration, and storage callbacks.
- `MEG`: complete derivative H/S assembly for one optimized basis function, with the
  same surrounding costs. It is not an isolated CUDA-kernel time.
- `EIG`: the generalized symmetric eigensolver. CPU rows use OpenBLAS/LAPACK;
  `GPU all` uses cuSOLVER.
- `wall`: process launch, initialization, input, all timed phases, and shutdown. For
  sub-millisecond matrix phases, wall time is dominated by MPI startup and I/O.

Speedup is always `CPU time / GPU time`; values above 1 favor the GPU.

## Physical workloads

| Code | System | Particles | Basis K | Source |
|---|---|---:|---:|---|
| `RG_0S` | Carbon, $^1S^e$ | 7 | 1000 | stored production basis |
| `RG_1P` | Lithium, $^2P^o$ | 4 | 250 | stored production basis |
| `RG_2D` | Lithium, $^2D^e$ | 4 | 100 | repository sample input |
| `RG_2P` | Carbon, $^3P^e$ | 7 | 600 | stored production basis |

These cases are representative calculations, not a cross-code microbenchmark. Different
symmetry projectors make equal K values unequal amounts of work.

## Fixed-basis H/S assembly: one V100 versus one CPU node

Jobs: GPU/CPU8 `7118`; CPU64 `7136`. All binaries and decks are byte-identical
between the two jobs. GPU MEE uses the `GPU ME` arm; wall uses `GPU all`.

| Code | CPU64 MEE (s) | V100 MEE (s) | MEE speedup | CPU64 EIG (s) | cuSOLVER EIG (s) | CPU64 wall (s) | GPU-all wall (s) | wall speedup |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `RG_0S` | 4.4978 | 0.9235 | **4.87x** | 0.1242 | 0.1004 | 9.336 | 3.510 | **2.66x** |
| `RG_1P` | 0.0074 | 0.0045 | **1.64x** | 0.0075 | 0.0715 | 4.477 | 2.138 | 2.09x* |
| `RG_2D` | 0.0014 | 0.0039 | **0.36x** | 0.0034 | 0.0711 | 4.578 | 2.624 | 1.74x* |
| `RG_2P` | 2.4797 | 0.9107 | **2.72x** | 0.0364 | 0.0948 | 7.053 | 3.459 | **2.04x** |

`*` The small-case wall ratios are launcher/I/O ratios, not accelerator ratios. Their MEE
rows give the meaningful crossover result.

The matrix path wins decisively for the symmetry-heavy Carbon cases, modestly for
`RG_1P` K=250, and loses for `RG_2D` K=100. cuSOLVER wins only at K=1000; at K=250,
K=100, and K=600 it is respectively 9.5x, 20.9x, and 2.6x slower than the CPU solve.
`ECG_GPU_EIG=1` should therefore remain opt-in.

## Eight-rank scaling context

| Code | CPU8 MEE (s) | V100 MEE (s) | MEE speedup | CPU8 EIG (s) | cuSOLVER EIG (s) | CPU8 wall (s) | GPU-all wall (s) |
|---|---:|---:|---:|---:|---:|---:|---:|
| `RG_0S` | 17.8684 | 0.9235 | **19.35x** | 0.1133 | 0.1004 | 22.516 | 3.510 |
| `RG_1P` | 0.0049 | 0.0045 | 1.09x | 0.0068 | 0.0715 | 4.350 | 2.138 |
| `RG_2D` | 0.0012 | 0.0039 | 0.31x | 0.0043 | 0.0711 | 4.489 | 2.624 |
| `RG_2P` | 12.8187 | 0.9107 | **14.08x** | 0.0340 | 0.0948 | 17.408 | 3.459 |

CPU64 improves full H/S throughput by 3.97x (`RG_0S`) and 5.17x (`RG_2P`) over
CPU8. It does not scale the sub-10-ms cases usefully because MPI overhead exceeds the
matrix work.

## Fixed-basis raw repetitions

Each list is execution order. Times are seconds.

### `RG_0S`, K=1000

| Arm | wall | MEE | EIG |
|---|---|---|---|
| CPU8 | 22.383, 22.491, 22.522, 22.516, 22.615 | 17.8536, 17.8684, 17.8982, 17.8272, 17.9093 | 0.1173, 0.1127, 0.1134, 0.1130, 0.1133 |
| CPU64 | 9.324, 9.336, 9.318, 9.351, 9.422 | 4.4958, 4.4938, 4.4978, 4.5025, 4.7249 | 0.1279, 0.1250, 0.1242, 0.1242, 0.1238 |
| GPU ME | 4.050, 3.481, 3.740, 3.687, 3.828 | 0.9273, 0.9136, 0.9235, 0.9152, 0.9312 | 0.2565, 0.2564, 0.2569, 0.2565, 0.2569 |
| GPU all | 3.501, 3.376, 3.559, 3.510, 3.550 | 0.9220, 0.9141, 0.9160, 0.9144, 0.9157 | 0.1255, 0.0985, 0.1004, 0.1015, 0.0990 |

`GPU ME` intentionally uses the CPU eigensolver inside the NVHPC/netlib binary; its EIG
column is not a competitive CPU baseline. `GPU all` replaces it with cuSOLVER.

### `RG_1P`, K=250

| Arm | wall | MEE | EIG |
|---|---|---|---|
| CPU8 | 4.607, 4.592, 4.350, 4.180, 4.302 | 0.0049, 0.0049, 0.0050, 0.0048, 0.0049 | 0.0081, 0.0082, 0.0067, 0.0068, 0.0068 |
| CPU64 | 4.464, 4.405, 4.477, 4.530, 6.162 | 0.0073, 0.0064, 0.0081, 0.0074, 0.0083 | 0.0075, 0.0077, 0.0073, 0.0077, 0.0074 |
| GPU ME | 2.552, 2.444, 2.110, 2.062, 2.082 | 0.0060, 0.0055, 0.0045, 0.0044, 0.0045 | 0.0064, 0.0066, 0.0062, 0.0063, 0.0062 |
| GPU all | 2.602, 2.548, 2.138, 2.106, 2.098 | 0.0053, 0.0051, 0.0042, 0.0044, 0.0047 | 0.0906, 0.0889, 0.0679, 0.0715, 0.0693 |

### `RG_2D`, K=100

| Arm | wall | MEE | EIG |
|---|---|---|---|
| CPU8 | 4.316, 4.545, 4.483, 4.539, 4.489 | 0.0012, 0.0012, 0.0012, 0.0012, 0.0012 | 0.0060, 0.0042, 0.0041, 0.0043, 0.0044 |
| CPU64 | 4.578, 4.601, 4.600, 4.504, 4.573 | 0.0014, 0.0015, 0.0014, 0.0014, 0.0014 | 0.0036, 0.0032, 0.0034, 0.0034, 0.0031 |
| GPU ME | 2.940, 2.499, 2.508, 2.425, 2.589 | 0.0039, 0.0040, 0.0038, 0.0038, 0.0039 | 0.0017, 0.0015, 0.0016, 0.0015, 0.0015 |
| GPU all | 2.661, 2.649, 2.624, 2.534, 2.583 | 0.0040, 0.0037, 0.0036, 0.0039, 0.0040 | 0.0919, 0.0708, 0.0691, 0.0711, 0.0711 |

### `RG_2P`, K=600

| Arm | wall | MEE | EIG |
|---|---|---|---|
| CPU8 | 17.408, 17.351, 17.431, 17.377, 17.426 | 12.9232, 12.8460, 12.7807, 12.8187, 12.8033 | 0.0341, 0.0340, 0.0342, 0.0336, 0.0335 |
| CPU64 | 7.097, 7.053, 7.071, 7.038, 6.884 | 2.4695, 2.4871, 2.4736, 2.4797, 2.4944 | 0.0364, 0.0364, 0.0365, 0.0360, 0.0363 |
| GPU ME | 3.459, 3.361, 3.506, 3.317, 3.612 | 0.9113, 0.9084, 0.9107, 0.9081, 0.9136 | 0.0591, 0.0608, 0.0593, 0.0598, 0.0585 |
| GPU all | 3.426, 3.459, 3.503, 3.435, 3.572 | 0.9112, 0.9078, 0.9091, 0.9087, 0.9093 | 0.0992, 0.0918, 0.0949, 0.0930, 0.0948 |

## Repeated derivative path

Job `7124`. Each cell starts from the same stored basis and performs one method-I
`OPT_CYCLE` over functions 1..8, at most 10 evaluations per function. All 24 cells did
exactly 89 energy builds and 8 gradient builds. cuSOLVER is disabled, so the backend is
the intended independent variable. Times below are milliseconds per build.

| Code | CPU8 MEE/build | V100 MEE/build | speedup | CPU8 MEG/build | V100 MEG/build | speedup | wall CPU8 / GPU (s) |
|---|---:|---:|---:|---:|---:|---:|---:|
| `RG_0S` | 235.180 | 13.079 | **17.98x** | 77.750 | 10.875 | **7.15x** | 26.324 / 4.072 |
| `RG_1P` | 0.112 | 0.157 | 0.71x | 0.125 | 0.125 | 1.00x | 4.343 / 2.398 |
| `RG_2D` | 0.045 | 0.135 | 0.33x | 0.125 | 0.250 | 0.50x | 4.281 / 2.470 |
| `RG_2P` | 185.472 | 13.382 | **13.86x** | 96.500 | 14.750 | **6.54x** | 21.766 / 3.887 |

The sub-millisecond `RG_1P` and `RG_2D` phase clocks are timer-resolution limited. They
show that these workloads are below crossover; their displayed ratios are not precision
measurements. A current-source CPU64 derivative campaign was not run.

### Derivative raw repetitions

| Code | Arm | wall (s) | MEE/build (ms) | MEG/build (ms) |
|---|---|---|---|---|
| `RG_0S` | CPU8 | 26.263, 26.488, 26.324 | 235.179775, 237.303371, 234.595506 | 79.375, 77.500, 77.750 |
| `RG_0S` | V100 | 4.177, 3.692, 4.072 | 13.078652, 13.078652, 13.067416 | 10.875, 10.875, 10.875 |
| `RG_1P` | CPU8 | 4.237, 4.343, 4.391 | 0.112360, 0.112360, 0.112360 | 0.125, 0.125, 0.125 |
| `RG_1P` | V100 | 2.386, 2.398, 2.475 | 0.157303, 0.157303, 0.146067 | 0.125, 0.125, 0.125 |
| `RG_2D` | CPU8 | 4.348, 4.231, 4.281 | 0.044944, 0.044944, 0.044944 | 0.125, 0.125, 0.125 |
| `RG_2D` | V100 | 2.478, 2.352, 2.470 | 0.134831, 0.134831, 0.134831 | 0.250, 0.250, 0.250 |
| `RG_2P` | CPU8 | 21.833, 21.766, 21.672 | 185.483146, 184.943820, 185.471910 | 96.500, 96.500, 97.000 |
| `RG_2P` | V100 | 3.888, 3.855, 3.887 | 13.382022, 13.359551, 13.404494 | 14.750, 14.750, 14.750 |

## Numerical agreement

Maximum absolute energy difference over all five fixed-basis repetitions:

| Code | CPU64 vs GPU matrix | CPU64 vs GPU matrix + cuSOLVER | default atomic GPU span |
|---|---:|---:|---:|
| `RG_0S` | 4.460e-10 | 7.585e-10 | 9.156e-10 |
| `RG_1P` | 1.776e-15 | 4.232e-12 | 0 |
| `RG_2D` | 0 | 7.629e-13 | 0 |
| `RG_2P` | 8.527e-14 | 1.631e-11 | 2.316e-11 |

Maximum CPU8/GPU endpoint difference in the repeated method-I derivative window:

| Code | max endpoint difference | three-run GPU endpoint span |
|---|---:|---:|
| `RG_0S` | 4.551e-10 | 6.224e-10 |
| `RG_1P` | 0 | 0 |
| `RG_2D` | 0 | 0 |
| `RG_2P` | 2.132e-14 | 0 |

Optimizer endpoints test self-consistency, not the gradient formula directly. A separate
analytic-versus-finite-difference gradient test remains open.

## Ordered H/S reduction

Job `7119` ran `ECG_DETERM=1` three times per code. All three full-precision energy
strings were identical for every code. Median MEE:

| Code | K | ordered V100 MEE (s) |
|---|---:|---:|
| `RG_0S` | 1000 | 0.8873 |
| `RG_1P` | 250 | 0.0052 |
| `RG_2D` | 100 | 0.0040 |
| `RG_2P` | 600 | 0.6902 |

`ECG_DETERM` applies only to energy-path H/S assembly. The derivative kernel retains
`atomicAdd` and optimization runs are not bit-reproducible. The ordered kernel needs
`16 * NumYHYTerms` bytes of dynamic shared memory per block, so it cannot launch when
the symmetry expansion exceeds the device limit (for example Oxygen, 40,320 terms).

## Validation inventory

| Job | Purpose | Valid cells | Result |
|---:|---|---:|---|
| `7114` | integrated four-code build/correctness gate | 28 | pass |
| `7118` | CPU8/V100 fixed-basis timings | 60 | pass |
| `7119` | ordered H/S reproducibility | 12 | pass |
| `7124` | repeated derivative path | 24 | pass |
| `7136` | full-node CPU64 timings | 20 | pass |
| | **Total** | **144** | **all pass** |

Job `7125` is not a numerical result. OpenMPI selected UCX/RDMA for a single-node
64-rank run and failed before ECGPACK started:

```text
UCX ERROR ibv_reg_mr(...) failed: Cannot allocate memory
Please set max locked memory ... current: 8192 kbytes
```

The corrected job `7136` raised the soft limit to 1 GiB and selected
`--mca pml ob1 --mca btl self,sm`. This is a benchmark-launcher fix, not an ECGPACK fix.

## Bounds on interpretation

- Only the four real-ECG energy codes are ported. `CG_0S` and the off-diagonal codes
  remain CPU-only.
- One V100 and one host were tested. Multi-GPU device placement is not performance-gated.
- `RG_1P` K=250 is near crossover; `RG_2D` K=100 is below it. Larger physical bases are
  required before assigning either code a general crossover K.
- The fixed-basis benchmark stops after H/S assembly and the method-G eigensolve. It does
  not time the CPU-only expectation-value/operator suite.
- Atomic GPU reductions make Carbon results stable to scientific precision but not bitwise
  reproducible. Use `ECG_DETERM=1` only within its stated energy-path/device-memory bounds.
