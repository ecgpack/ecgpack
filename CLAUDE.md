# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

ECGPACK is a collection of closely related parallel Fortran codes for variational calculations of quantum few-body systems (few-electron atoms and molecules) using all-particle explicitly correlated Gaussian (ECG) basis functions. Developed by the research group of Sergiy Bubin (Physics Department, Nazarbayev University). Parallelism is via MPI. The public repository is <https://github.com/ecgpack/ecgpack> and its default branch is `main`. See `README.md` and directory `doc/` for documentation that includes physics conventions and references.

## Repository layout

Each code lives in its own top-level directory with an identical structure: a `Makefile`, a `src/` subdirectory, and a `.vscode/` subdirectory (nearly identical, user-independent JSON configs for building/debugging the code in VSCode). Most, but not all, codes also include a `sample_input/` subdirectory of worked examples (see below). The codes fall into several groups (some of these groups may overlap):

- **Energy codes (diagonal matrix element codes; basis generation codes)** — `RG_0S`, `RG_1P`, `RG_2D`, `RG_2P`, `CG_0S`. These are used for generating ECG basis sets for states of different angular-momentum/parity. They can compute energies, expectation values, particle distributions, and save wave functions.
- **Off-diagonal matrix element codes** — `RG_0S-1P`, `RG_0S-2D`, `RG_0S-2P`, `RG_1P-1P`, `RG_1P-2D`, `RG_1P-2P`, `RG_2D-2D`, `RG_2P-2D`, `RG_2P-2P`. These codes are used to compute off-diagonal matrix elements between two states of different angular-momentum/parity (and which are thus expanded using different basis sets) for such operators as transition dipole moment, spin–orbit interaction, and non-contact spin–spin interaction.
- **Transition dipole codes** — `RG_0S-1P`, `RG_1P-2D`, `RG_1P-2P`. These codes are used to compute off-diagonal matrix elements of the transition dipole moment operator in the length and velocity gauges.
- **Fine structure coupling codes** — `RG_0S-2D`, `RG_0S-2P`, `RG_1P-1P`, `RG_2D-2D`, `RG_2P-2D`, `RG_2P-2P`. These codes are used to compute off-diagonal matrix elements of the spin–orbit interaction and non-contact spin–spin interaction.
- **Real Gaussian codes (Real ECG codes)** — `RG_0S`, `RG_1P`, `RG_2D`, `RG_2P`, `RG_0S-1P`, `RG_0S-2D`, `RG_0S-2P`, `RG_1P-1P`, `RG_1P-2D`, `RG_1P-2P`, `RG_2D-2D`, `RG_2P-2D`, `RG_2P-2P`. These deal with real Gaussian basis sets (real ECGs).
- **Complex Gaussian codes (Complex ECG codes)** — `CG_0S`. These deal with complex ECG basis sets (complex ECGs). Currently there is only one basis ($L=0$), and it is a work in progress that lags behind the real ECG codes. It may not even be compilable.

Most codes ship with a `sample_input/` subdirectory of worked examples — the four real-ECG energy codes (`RG_0S`, `RG_1P`, `RG_2D`, `RG_2P`) plus all nine off-diagonal codes (`RG_0S-1P`, `RG_0S-2D`, `RG_0S-2P`, `RG_1P-1P`, `RG_1P-2D`, `RG_1P-2P`, `RG_2D-2D`, `RG_2P-2D`, `RG_2P-2P`). Only `CG_0S` currently has no sample inputs. See the Running section below.

`RG_0S` is the most complete reference implementation; other codes share much of its source and Makefile. `CG_0S/list_of_things_to_do.txt` holds older development notes — check each item against the current sources before acting on it, since some (e.g. precomputing scaled charge products) are already done. Non-code directories: `doc/` (documentation: `compilation_and_execution.md`, `input_file_format.md`, `theoretical_background.md`, `use_of_visual_studio_code.md`, plus developer-facing design guides in `doc/devnotes/` and logos in `doc/img/`), `utilities/` (Python/Bash helpers for input-file generation and merging, nuclear masses, job monitoring, and a Young-tableaux spreadsheet), and `bin/` + `jobs/` (created by the user, not in git — both are gitignored, as are the per-code `debug/`, `debug-omp/`, `release/`, `release-omp/` trees).

The root also holds `README.md` (the repository manual), `AGENTS.md` (general guidance for AI coding agents — **keep it and this file consistent when either changes**), `AUTHORS.md` (contributors), `CITATION.cff` (citation metadata), `LICENSE.md` (BSD 3-Clause), `THIRD-PARTY-NOTICES.md` (notices for the bundled qrupdate-ng, BLAS, LAPACK, PORT, and SLATEC-derived sources — keep it in sync when bundled third-party code changes), `build.bash`, and `.code-workspace` (a VSCode multi-folder workspace grouping all the code directories).

## Building

The canonical entry point is `build.bash` in the root directory — run it from that directory, with no arguments, to print full usage (it then exits 1). It loops over toolchains/configs/codes/nparticles/precisions/linalg/openmp choices, builds via each code's Makefile, and stores binaries in `bin/[<machine>/]<toolchain>/<config>/<CODE>_N<nparticles>_P<precision>_<linalg>[_omp]`. Notes on that path and name: the `<machine>/` segment is inserted only for the cluster machines `puma`, `ocelote`, and `elgato` (it is empty for `linux-generic`, `ubuntu-generic`, `irgetas`, `shabyt`, `muon`); the `<linalg>` suffix is always present (e.g. `_netlib`, `_mkl`, `_openblas`); OpenMP builds add a trailing `_omp` but land in the same `<config>` directory as the serial ones.

Only `nparticles` is required. Defaults: `machine=linux-generic`, `toolchain=systemdefault`, `config=debug,release`, `code=` all codes **except** `CG_0S`, `precision=8`, `linalg=netlib`, `openmp=no`. Every argument except `machine` accepts a comma-separated list. The `linalg` argument (see below) takes one of `netlib` (default), `mkl`, `lblas`, `openblas`, `aocl`; for `precision=10`/`16` only `netlib` is built (other values are skipped), and `openmp=yes` is skipped for any code other than `RG_0S`/`RG_1P`/`RG_2D`/`RG_2P`.

Supported machines: `linux-generic` (default), `ubuntu-generic`, `irgetas`, `shabyt`, `muon`, `puma`, `ocelote`, `elgato`. Toolchains are per machine — the generic/NU machines take `systemdefault`, `foss-2023b`/`2024a`/`2025a`/`2025b`, `intel-2023b`/`2024a`/`2025a`/`2025b`, `nvhpc-25.3`/`25.9`/`26.5`; `puma` takes `systemdefault`, `ohpc`, `intel-2024.0.0`; `ocelote`/`elgato` take `systemdefault`, `ohpc`, `intel-2020.4`. Non-`systemdefault` toolchains load Easybuild/`module` environments. Machine profiles are hardcoded in both `build.bash` and the Makefiles (adding one means editing both), but toolchain names and module-loading logic live only in `build.bash` — the Makefiles pick flags from `MACHINE`, `COMPILER`, and `PREC`.

```bash
./build.bash machine=ubuntu-generic toolchain=systemdefault config=release code=RG_0S nparticles=4 precision=8 linalg=netlib openmp=no
```

To build a single code directly, invoke its Makefile (this is what `build.bash` calls under the hood):

```bash
cd RG_0S
make release COMPILER=gfortran MACHINE=ubuntu-generic PREC=8 LINALG=openblas EXEFILE=ecg
make release COMPILER=gfortran MACHINE=ubuntu-generic PREC=8 LINALG=openblas OPENMP=yes EXEFILE=ecg
make debug   COMPILER=gfortran MACHINE=ubuntu-generic PREC=8 LINALG=netlib   EXEFILE=ecg
make clean   # also: cleaner, cleanest, cleanrelease, cleandebug
```

A direct Makefile build leaves the executable in its own build tree under the `EXEFILE` name — `release/ecg`, or `release-omp/ecg` with `OPENMP=yes`. `COMPILER`, `PREC`, and `LINALG` all share one object directory, so `make clean` before switching any of them; `build.bash` sidesteps this by running `make clean` for every build combination.

Key build parameters:

- `CONFIG` (`release`/`debug`) — set by the `make release`/`make debug` target. `release` optimizes (`-O3 -march=native` for gfortran, `-O3 -xhost -fp-model precise` for ifort/ifx, `-O3 -tp=native` for nvfortran); `debug` enables bounds/uninit/FPE checks (`-O0 -g -fcheck=all -ffpe-trap=...` and equivalents). Object and `.mod` files go in `release/` or `debug/`.
- `PREC` — real `kind`: `8` (double/fp64), `10` (extended/fp80), `16` (quadruple). Selects which `src/wp_def_<PREC>.f90` is compiled. The Makefiles mark as supported: `8`/`10`/`16` for gfortran, `8`/`16` for ifort and ifx, `8` only for nvfortran. Compile support is not run support — gfortran + `PREC=16` currently runs **serial only** (an OpenMPI limitation); see `doc/compilation_and_execution.md` for the full matrix.
- `LINALG` — selects which BLAS/LAPACK implementation to link against (default `netlib`). `netlib` compiles the bundled, lightly modified reference `src/BLAS.f`/`src/LAPACK.f` and adds no extra link flags; `mkl` (Intel MKL — compiler-dependent `-lmkl_*` flags), `lblas` (`-llapack -lblas`), `openblas` (`-lopenblas`), and `aocl` (AMD AOCL `-lflame -lblis …`) instead link an external optimized library and skip the bundled sources. Only `PREC=8` honors the optimized choices; for `PREC=10`/`16` only `LINALG=netlib` is supported (any other value leaves the build unsupported). In the off-diagonal codes `LINALG` is accepted but a no-op (they link no BLAS/LAPACK). The bundled `src/BLAS.f`/`src/LAPACK.f` and the `BARE_OBJS_LPKBLS` object list are compiled only when `LINALG=netlib`.
- `OPENMP` (`no`/`yes`, default `no`) — supported by the four real-ECG energy codes; adds the compiler's OpenMP flag (`-fopenmp` gfortran, `-qopenmp` ifort/ifx, `-mp` nvfortran) and isolates objects in `debug-omp/` or `release-omp/`. The actual `!$OMP` directives live in those codes' bundled `src/BLAS.f`, gated by an `OMP_MIN_WORK` size threshold, and so pay off mainly in `LINALG=netlib` builds through qrlinalg; an external BLAS brings its own threading instead. At runtime set `OMP_NUM_THREADS` explicitly and budget roughly (MPI ranks x threads per rank) cores — each rank spawns its own team, and nested threading from an external BLAS can oversubscribe a node.
- `COMPILER` (`gfortran`→`mpif90`, `ifort`→`mpiifort`, `ifx`→`mpiifx`, `nvfortran`→`mpif90`) and `MACHINE` together select the compiler flags (see the machine/toolchain notes above).

Note: the **number of particles is compiled in**, not a runtime argument. `build.bash` does an in-place `sed` on `src/wp_def_<PREC>.f90` to set `Glob_AllowedNumOfParticles` and `MPI_DPREC`, builds, then restores the original from `wp_def_temporary.f90`; it accepts `nparticles` from 2 to 19. A direct Makefile build uses whatever values are already in that file. A binary built for N particles rejects input files with a different `PARTICLES` count. Do not run two builds of the same code directory concurrently — they share the backup file and build trees — and after an interrupted build, check that `src/wp_def_<PREC>.f90` was restored.

## Running

Each binary is an MPI program. The energy codes read a single input/output file named **`inout.txt`** from the current working directory (default `Glob_DataFileName` in `globvars.f90`); the file's first (optional) line is `BASIS_TYPE <BASIS>`, the second (required) line is `PARTICLES <n>`. The off-diagonal matrix-element codes instead read two wave-function files (one per state) and do not use `inout.txt`. The transition dipole codes (`RG_0S-1P`, `RG_1P-2D`, `RG_1P-2P`), as well as `RG_0S-2D`, read just the two wave-function files directly, with no extra input file; the remaining fine-structure coupling codes (`RG_0S-2P`, `RG_1P-1P`, `RG_2D-2D`, `RG_2P-2D`, `RG_2P-2P`) additionally require a small `inp.txt` whose positional argument(s) select the transition/coupling type. The wave-function file names vary by code: most use `wf_state0.txt`/`wf_state1.txt`, but `RG_1P-2D` and `RG_1P-2P` use `wf_state1.txt`/`wf_state2.txt`. Run from a job directory containing the required input file(s):

```bash
mpirun -np <N> /path/to/bin/[<machine>/]<toolchain>/<config>/RG_0S_N4_P8_netlib
```

When reproducing a sample case, copy it into `jobs/` (or another scratch directory) and run it there with an absolute path to the binary — the codes rewrite `inout.txt` in place, so running inside `sample_input/` overwrites tracked examples.

There is no automated test suite; validation is done by running physical test cases and comparing computed energies/expectation values against published reference values. The codes that ship sample inputs (see Repository layout) provide ready-to-run cases under `<CODE>/sample_input/`, documented by `README.md` files giving the description, expected runtime, and the reference energy/values to reproduce. Where those READMEs sit depends on the code: the energy and transition dipole codes put one in each case directory, whereas the fine-structure codes put one in each of the case's three subdirectories (`initial_state_*/`, `final_state_*/`, `matelem/`) and none at the case level. `RG_0S-1P` and `RG_1P-2D` additionally have an overview `README.md` at the top of `sample_input/`. For the energy codes each case is a single directory with an `inout.txt`; recurring kinds are `basis_generation_*` (build an ECG basis from scratch up to a target size), `expected_values_*` (compute expectation values from a saved basis), `densities_*` (compute densities and pair correlation functions in coordinate and momentum space via the `DENSITIES`/`MOMT_DENS` steps, shipping Python grid builders — `create_grid.py`, or `create_radial_grid.py` + `create_cyl_grid.py` — a sample grid `.dat`, gnuplot plotting scripts `create_plots_{coor,mom}_{dens,cf}.gp`, and sample output/plots in `output_cyl/` and `cyl_plots/`), and `store_wavefunction_*` (saves both the linear and nonlinear variational parameters of the wave function into a file). The off-diagonal codes' examples instead bundle subdirectories for the two states plus the matrix-element run: the transition dipole cases (`RG_0S-1P`, `RG_1P-2D`, `RG_1P-2P`) are named `transition_dipole_moment_*`, while the fine-structure cases (`RG_0S-2D`, `RG_0S-2P`, `RG_1P-1P`, `RG_2D-2D`, `RG_2P-2D`, `RG_2P-2P`) provide `initial_state_*/`, `final_state_*/`, and `matelem/` subdirectories (the `matelem/` directory holds the two wave-function files and, for every fine-structure code except `RG_0S-2D`, the `inp.txt` that selects the coupling type).

## Source architecture

In the four real-ECG energy codes (the fullest source set), the module compile/dependency order (see the Makefile) is:

`wp_def_<PREC>` → `globvars` → `misc`, `linalg`, `spin` → `matelem` → `matform` → `qrupdate` → `qrlinalg` → `workproc` → `main`

The other codes carry a subset of these files:

- **`CG_0S`** — same as above minus `spin.f90`, `qrlinalg.f90`, and `qrupdate/`. It also lacks the `'Q'` eigenvalue method and the `MOMT_DENS`/`SAVE_HSWF` steps, so do not assume parity with the `RG_*` energy codes.
- **Off-diagonal codes** — only `wp_def_<PREC>`, `globvars`, `matelem`, (`spin` where spin coupling is needed), `workproc`, `main`. No `misc`/`linalg`/`matform`/`qrupdate`, and no bundled BLAS/LAPACK. `spin.f90` is present in `RG_0S-2D`, `RG_0S-2P`, `RG_1P-1P`, `RG_2D-2D`, `RG_2P-2D`, `RG_2P-2P` and absent in the transition-dipole codes `RG_0S-1P`, `RG_1P-2D`, `RG_1P-2P`.

Key files:

- **`wp_def_<PREC>.f90`** — defines the working-precision kind (`wp`), `Glob_AllowedNumOfParticles`, and the MPI real type. Edited at build time by `build.bash` (see above).
- **`globvars.f90`** — all global state (the `Glob_*` variables: masses, charges, basis, matrices) and physical/numeric constants.
- **`linalg.f90`** — linear-algebra wrappers over BLAS/LAPACK, including serial/MPI routing and cache-blocked LDL factorization paths. `linalg_setparam(n)` measures and selects those paths once per problem dimension; it must stay a **collective call on all MPI ranks** and be called again whenever the dimension changes. `BLAS.f`/`LAPACK.f` are bundled, lightly modified netlib reference sources used only when `LINALG=netlib` (the qrupdate subdirectory no longer carries its own duplicate copies — it shares `src/BLAS.f`/`src/LAPACK.f`). `dmng.f` (lightly modified TOMS nonlinear minimizer) and `X1MACH.f90` (machine constants) support the optimizer.
- **`spin.f90`** — spin algebra and permutational-symmetry projection.
- **`matelem.f90`** — matrix elements between individual basis functions; **`matform.f90`** assembles the full Hamiltonian (H) and overlap (S) matrices. Energy-code entry points are symmetry-qualified: `MatrixElementsHS_<BASIS>` (H/S, optionally with nonlinear-parameter derivatives), `MatrixElementsAll_<BASIS>` (full operator/expectation-value set), and `OverlapMatrixElement_<BASIS>` / `NormalizedOverlapMatElem_<BASIS>` — e.g. `MatrixElementsHS_RG_0S`, `MatrixElementsAll_CG_0S`. Off-diagonal codes use their own operator-specific routines; read the local `matelem.f90` rather than assuming a name. This is performance-sensitive code: preserve the preassembled matrix expressions, the no-gradient fast paths, and the use of the precomputed `Glob_ScaledPseudoChargeMatrix` (built once in `ProgramDataInit`) instead of recomputing scaled charge products in hot loops.
- **`qrlinalg.f90` and `qrupdate/`** — QR factorization/update support for generalized-eigenvalue method `'Q'` in the four real-ECG energy codes. The contract: basis arrays stay in canonical order, and the **lower** triangles of `Glob_H`/`Glob_S` (including diagonals) are authoritative, holding normalized, unshifted physical matrix elements — the upper triangles are scratch/unspecified. Read `doc/devnotes/Q_method_design.md` before touching Q matrix assembly, updates, or solver state; it also documents the porting procedure.
- **`workproc.f90`** — the bulk of the program (~13,000–14,000 lines in the real-ECG energy codes, ~9,000 in `CG_0S`, ~1,300–2,200 in the off-diagonal codes): `ReadIOFile`/`SaveResults` I/O, basis enlargement, optimization cycles, the generalized symmetric eigenvalue solvers (methods `'G'`, `'I'`, and `'Q'` in the real-ECG energy codes), expectation values, densities, and swap-file handling.
- **`main.f90`** — initializes MPI, seeds RNGs, then drives a sequence of **BBOP** (Basis Building and Optimization Program) steps read from the input file. Each step is a `case` in the main `select`: `BASIS_ENL`, `OPT_CYCLE`, `FULL_OPT1`, `ELIM_LCFN`, `ELIM_LND1`, `SEPR_LND1`, `SEPR_FLCF`, `EXPC_VALS`, `DENSITIES`, `MOMT_DENS`, `SAVE_FILE`, `SAVE_HSWF`. Adding a calculation mode means adding a case here plus the corresponding routine in `workproc.f90`.

When editing matrix-element or matform code, changes usually must be mirrored across the analogous `RG_*` directories, since the codes are intentionally near-duplicates specialized to different symmetries. The same applies to shared build logic, file formats, and solver changes. Prefer small, targeted edits in the existing style; preserve scientific behavior (and third-party license banners in bundled sources) unless a change is explicitly requested, and validate numerical changes against a sample case when practical.

## Editing and validation

Before editing, work out which code directories are affected and whether the same pattern exists in sibling `RG_*` directories; check the relevant `doc/` file when touching input formats, build behavior, or physics conventions. After editing, run a focused `make debug` (or `make release`) for the affected code, run a small sample case when numerical behavior could change, and confirm nothing leaked into `debug*/`, `release*/`, `bin/`, `jobs/`, or `src/wp_def_temporary.f90`. Keep this file and `AGENTS.md` aligned when changing agent-facing guidance.

## Further documentation

`doc/compilation_and_execution.md` (build/run details), `doc/input_file_format.md` (the `inout.txt` format and every BBOP step's arguments), `doc/theoretical_background.md` (physics and notation), `doc/use_of_visual_studio_code.md` (VSCode setup), `doc/devnotes/Q_method_design.md` (the `'Q'` method's canonical H/S layout and porting procedure).
