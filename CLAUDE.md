# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

ECG is a collection of closely related parallel Fortran codes for variational calculations of quantum few-body systems (few-electron atoms and molecules) using all-particle explicitly correlated Gaussian (ECG) basis functions. Developed by the research group of Sergiy Bubin (Physics Department, Nazarbayev University). Parallelism is via MPI. See `README.md` and directory `doc/` for documentation that includes physics conventions and references.

## Repository layout

Each code lives in its own top-level directory with an identical structure: a `Makefile`, a `src/` subdirectory, and a `.vscode/` subdirectory (nearly identical, user-independent JSON configs for building/debugging the code in VSCode). Most, but not all, codes also include a `sample_input/` subdirectory of worked examples (see below). The codes fall into several groups (some of these groups may overlap):

- **Energy codes (diagonal matrix element codes; basis generation codes)** — `RG_0S`, `RG_1P`, `RG_2D`, `RG_2P`, `CG_0S`. These are used for generating ECG basis sets for states of different angular-momentum/parity. They can compute energies, expectation values, particle distributions, and save wave functions.
- **Off-diagonal matrix element codes** — `RG_0S-1P`, `RG_0S-2D`, `RG_0S-2P`, `RG_1P-1P`, `RG_1P-2D`, `RG_1P-2P`, `RG_2D-2D`, `RG_2P-2D`, `RG_2P-2P`. These codes are used to compute off-diagonal matrix elements between two states of different angular-momentum/parity (and which are thus expanded using different basis sets) for such operators as transition dipole moment, spin–orbit interaction, and non-contact spin–spin interaction.
- **Transition dipole codes** — `RG_0S-1P`, `RG_1P-2D`, `RG_1P-2P`. These codes are used to compute off-diagonal matrix elements of the transition dipole moment operator in the length and velocity gauges.
- **Fine structure coupling codes** — `RG_0S-2D`, `RG_0S-2P`, `RG_1P-1P`, `RG_2D-2D`, `RG_2P-2D`, `RG_2P-2P`. These codes are used to compute off-diagonal matrix elements of the spin–orbit interaction and non-contact spin–spin interaction.
- **Real Gaussian codes (Real ECG codes)** — `RG_0S`, `RG_1P`, `RG_2D`, `RG_2P`, `RG_0S-1P`, `RG_0S-2D`, `RG_0S-2P`, `RG_1P-1P`, `RG_1P-2D`, `RG_1P-2P`, `RG_2D-2D`, `RG_2P-2D`, `RG_2P-2P`. These deal with real Gaussian basis sets (real ECGs).
- **Complex Gaussian codes (Complex ECG codes)** — `CG_0S`. These deal with complex ECG basis sets (complex ECGs). Currently there is only one basis ($L=0$), and it is a work in progress that lags behind the real ECG codes. It may not even be compilable.

Most codes ship with a `sample_input/` subdirectory of worked examples — the four real-ECG energy codes (`RG_0S`, `RG_1P`, `RG_2D`, `RG_2P`) plus all nine off-diagonal codes (`RG_0S-1P`, `RG_0S-2D`, `RG_0S-2P`, `RG_1P-1P`, `RG_1P-2D`, `RG_1P-2P`, `RG_2D-2D`, `RG_2P-2D`, `RG_2P-2P`). Only `CG_0S` currently has no sample inputs. See the Running section below.

`RG_0S` is the most complete reference implementation; other codes share much of its source and Makefile. Non-code directories: `doc/` (documentation), `utilities/`, and `bin/` + `jobs/` (created by the user, not in git). The root also holds `README.md` (the repository manual), `AUTHORS.md` (list of contributors), and `.code-workspace` (a VSCode multi-folder workspace grouping all the code directories).

## Building

The canonical entry point is `build.bash` in the root directory — run it with no arguments for full usage. It loops over toolchains/configs/codes/precisions/linalg choices, builds via each code's Makefile, and stores binaries in `bin/<toolchain>/<config>/<CODE>_N<nparticles>_P<precision>_<linalg>` (the `<linalg>` suffix is always present — e.g. `_netlib`, `_mkl`, `_openblas`). The `nparticles` argument is required. The `linalg` argument (see below) takes one of `netlib` (default), `mkl`, `lblas`, `openblas`, `aocl`; for `precision=10`/`16` only `netlib` is built (other values are skipped).

```bash
./build.bash machine=ubuntu-generic toolchain=systemdefault config=release code=RG_0S nparticles=4 precision=8
```

To build a single code directly, invoke its Makefile (this is what `build.bash` calls under the hood):

```bash
cd RG_0S
make release COMPILER=gfortran MACHINE=ubuntu-generic PREC=8 LINALG=openblas EXEFILE=ecg
make debug   COMPILER=gfortran MACHINE=ubuntu-generic PREC=8 LINALG=netlib   EXEFILE=ecg
make clean   # also: cleaner, cleanest, cleanrelease, cleandebug
```

Key build parameters:

- `CONFIG` (`release`/`debug`) — set by the `make release`/`make debug` target. `release` uses `-O3 -march=native`; `debug` enables bounds/uninit/FPE checks. Object and `.mod` files go in `release/` or `debug/`.
- `PREC` — real `kind`: `8` (double/fp64), `10` (extended/fp80, GNU only), `16` (quadruple). Selects which `src/wp_def_<PREC>.f90` is compiled.
- `LINALG` — selects which BLAS/LAPACK implementation to link against (default `netlib`). `netlib` compiles the bundled, lightly modified reference `src/BLAS.f`/`src/LAPACK.f` and adds no extra link flags; `mkl` (Intel MKL — compiler-dependent `-lmkl_*` flags), `lblas` (`-llapack -lblas`), `openblas` (`-lopenblas`), and `aocl` (AMD AOCL `-lflame -lblis …`) instead link an external optimized library and skip the bundled sources. Only `PREC=8` honors the optimized choices; for `PREC=10`/`16` only `LINALG=netlib` is supported (any other value leaves the build unsupported). In the off-diagonal codes `LINALG` is accepted but a no-op (they link no BLAS/LAPACK). The bundled `src/BLAS.f`/`src/LAPACK.f` and the `BARE_OBJS_LPKBLS` object list are compiled only when `LINALG=netlib`.
- `COMPILER` (`gfortran`→`mpif90`, `ifort`→`mpiifort`, `ifx`→`mpiifx`, `nvfortran`→`mpif90`) and `MACHINE` select compiler flags. Supported machines are hardcoded in both `build.bash` and the Makefiles; adding a machine means editing both.

Note: the **number of particles is compiled in**, not a runtime argument. `build.bash` does an in-place `sed` on `src/wp_def_<PREC>.f90` to set `Glob_AllowedNumOfParticles`, builds, then restores the original from `wp_def_temporary.f90`. A binary built for N particles rejects input files with a different particle count.

## Running

Each binary is an MPI program. The energy codes read a single input/output file named **`inout.txt`** from the current working directory (default `Glob_DataFileName` in `globvars.f90`); the file's first (optional) line is `BASIS_TYPE <BASIS>`, the second (required) line is `PARTICLES <n>`. The off-diagonal matrix-element codes instead read two wave-function files (one per state) and do not use `inout.txt`. The transition dipole codes (`RG_0S-1P`, `RG_1P-2D`, `RG_1P-2P`), as well as `RG_0S-2D`, read just the two wave-function files directly, with no extra input file; the remaining fine-structure coupling codes (`RG_0S-2P`, `RG_1P-1P`, `RG_2D-2D`, `RG_2P-2D`, `RG_2P-2P`) additionally require a small `inp.txt` whose positional argument(s) select the transition/coupling type. The wave-function file names vary by code: most use `wf_state0.txt`/`wf_state1.txt`, but `RG_1P-2D` and `RG_1P-2P` use `wf_state1.txt`/`wf_state2.txt`. Run from a job directory containing the required input file(s):

```bash
mpirun -np <N> /path/to/bin/<toolchain>/<config>/RG_0S_N4_P8_netlib
```

There is no automated test suite; validation is done by running physical test cases and comparing computed energies/expectation values against published reference values. The codes that ship sample inputs (see Repository layout) provide ready-to-run cases under `<CODE>/sample_input/`, documented by `README.md` files giving the description, expected runtime, and the reference energy/values to reproduce. Where those READMEs sit depends on the code: the energy and transition dipole codes put one in each case directory, whereas the fine-structure codes put one in each of the case's three subdirectories (`initial_state_*/`, `final_state_*/`, `matelem/`) and none at the case level. `RG_0S-1P` and `RG_1P-2D` additionally have an overview `README.md` at the top of `sample_input/`. For the energy codes each case is a single directory with an `inout.txt`; recurring kinds are `basis_generation_*` (build an ECG basis from scratch up to a target size), `expected_values_*` (compute expectation values from a saved basis), `densities_*` (compute densities and pair correlation functions, with grid-builder Bash scripts, gnuplot plotting scripts, and sample output plots), and `store_wavefunction_*` (saves both the linear and nonlinear variational parameters of the wave function into a file). The off-diagonal codes' examples instead bundle subdirectories for the two states plus the matrix-element run: the transition dipole cases (`RG_0S-1P`, `RG_1P-2D`, `RG_1P-2P`) are named `transition_dipole_moment_*`, while the fine-structure cases (`RG_0S-2D`, `RG_0S-2P`, `RG_1P-1P`, `RG_2D-2D`, `RG_2P-2D`, `RG_2P-2P`) provide `initial_state_*/`, `final_state_*/`, and `matelem/` subdirectories (the `matelem/` directory holds the two wave-function files and, for every fine-structure code except `RG_0S-2D`, the `inp.txt` that selects the coupling type).

## Source architecture

Within each `src/`, the module compile/dependency order (see the Makefile) is:

`wp_def_<PREC>` → `globvars` → `misc`, `linalg`, `spin` → `matelem` → `matform` → `workproc` → `main`

- **`wp_def_<PREC>.f90`** — defines the working-precision kind (`wp`), `Glob_AllowedNumOfParticles`, and the MPI real type. Edited at build time by `build.bash` (see above).
- **`globvars.f90`** — all global state (the `Glob_*` variables: masses, charges, basis, matrices) and physical/numeric constants.
- **`linalg.f90`** — linear-algebra wrappers over BLAS/LAPACK. `BLAS.f`/`LAPACK.f` are bundled, lightly modified netlib reference sources used only when `LINALG=netlib`. `dmng.f` (lightly modified TOMS nonlinear minimizer) and `X1MACH.f90` (machine constants) support the optimizer.
- **`spin.f90`** — spin algebra and permutational-symmetry projection.
- **`matelem.f90`** — matrix elements between individual basis functions; **`matform.f90`** assembles the full Hamiltonian (H) and overlap (S) matrices.
- **`workproc.f90`** — the bulk of the program (often >8000 lines): `ReadIOFile`/`SaveResults` I/O, basis enlargement, optimization cycles, the generalized symmetric eigenvalue solvers (methods `'G'` and `'I'`), expectation values, densities, and swap-file handling.
- **`main.f90`** — initializes MPI, seeds RNGs, then drives a sequence of **BBOP** (Basis Building and Optimization Program) steps read from the input file. Each step is a `case` in the main `select`: `BASIS_ENL`, `OPT_CYCLE`, `FULL_OPT1`, `ELIM_LCFN`, `ELIM_LND1`, `SEPR_LND1`, `SEPR_FLCF`, `EXPC_VALS`, `DENSITIES`, `MOMT_DENS`, `SAVE_FILE`, `SAVE_HSWF`. Adding a calculation mode means adding a case here plus the corresponding routine in `workproc.f90`.

When editing matrix-element or matform code, changes usually must be mirrored across the analogous `RG_*` directories, since the codes are intentionally near-duplicates specialized to different symmetries.
