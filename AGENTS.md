# AGENTS.md

Guidance for AI coding agents working in this repository.

## Project Overview

ECGPACK is a collection of closely related parallel Fortran codes for high-accuracy variational calculations of quantum few-body systems, including few-electron atoms, molecules, ions, and systems with exotic particles. The codes use all-particle explicitly correlated Gaussian (ECG) basis functions and MPI parallelism. The public repository is <https://github.com/ecgpack/ecgpack>, and its default branch is `main`.

The project is developed by the research group of Sergiy Bubin, Physics Department, Nazarbayev University. For scientific background, notation, and references, see `README.md` and the documents in `doc/`.

## Repository Layout

Each main code lives in a top-level directory with a similar structure:

- `Makefile`
- `src/`
- `.vscode/`: per-code build tasks, debugger launch configurations, and editor settings; see `doc/use_of_visual_studio_code.md`
- usually `sample_input/`

The main code groups overlap:

- Energy and wave-function codes: `RG_0S`, `RG_1P`, `RG_2D`, `RG_2P`, `CG_0S`. These generate ECG basis sets for states of different angular momentum and parity, and can compute energies, expectation values, particle distributions, and saved wave functions.
- Off-diagonal matrix-element codes: `RG_0S-1P`, `RG_0S-2D`, `RG_0S-2P`, `RG_1P-1P`, `RG_1P-2D`, `RG_1P-2P`, `RG_2D-2D`, `RG_2P-2D`, `RG_2P-2P`. These work with pairs of states expanded in different symmetry-adapted bases.
- Transition dipole codes: `RG_0S-1P`, `RG_1P-2D`, `RG_1P-2P`. These evaluate transition electric-dipole moments in the length and velocity gauges.
- Fine-structure coupling codes: `RG_0S-2D`, `RG_0S-2P`, `RG_1P-1P`, `RG_2D-2D`, `RG_2P-2D`, `RG_2P-2P`. These evaluate spin-orbit and non-contact spin-spin matrix elements.
- Real ECG codes: all `RG_*` directories.
- Complex ECG codes: `CG_0S`. This code currently supports only the complex $L=0$ basis and is a work in progress that may lag behind the real-ECG codes.

`RG_0S` is the most complete reference implementation. Other codes intentionally share much of its source layout and Makefile structure.

`CG_0S/list_of_things_to_do.txt` contains older development notes. Check them against the current sources before acting; some listed work, such as precomputing scaled charge products, is already implemented.

Non-code directories:

- `doc/`: project documentation (`compilation_and_execution.md`, `input_file_format.md`, `theoretical_background.md`, `use_of_visual_studio_code.md`); `doc/devnotes/` contains developer-facing design and maintenance guides, and `doc/img/` holds the project logos
- `utilities/`: helper scripts — input-file generation and merging (`input_file_manipulation/`), nuclear masses (`nuclear_masses/`), job monitoring (`job_monitoring/`), and a Young-tableaux spreadsheet (`permutational_symmetry/`)
- `bin/`: generated/user-created binaries, not normally committed
- `jobs/`: user-created calculation work directories, not normally committed

`bin/`, `jobs/`, and the per-code `debug/`, `debug-omp/`, `release/`, `release-omp/` trees are gitignored, as is the temporary `*/src/wp_def_temporary.f90` that `build.bash` creates during a build.

Root files of interest:

- `README.md`: main project manual
- `AUTHORS.md`: contributor list
- `CITATION.cff`: machine-readable citation metadata
- `LICENSE.md`: BSD 3-Clause license for ECGPACK
- `THIRD-PARTY-NOTICES.md`: notices for bundled qrupdate-ng, BLAS, LAPACK, PORT, and SLATEC-derived sources
- `build.bash`: batch build driver
- `.code-workspace`: VS Code multi-folder workspace
- `CLAUDE.md`: Claude-specific project context; keep it consistent with this file when relevant

## General Agent Rules

- Prefer small, targeted changes that follow the existing code style.
- Do not rewrite broad portions of the duplicated Fortran code unless the task explicitly requires it.
- Treat the `RG_*` directories as related implementations. When changing shared algorithms, matrix elements, build logic, or file formats, check whether analogous changes are needed in sibling directories.
- Do not commit generated build products from `debug/`, `debug-omp/`, `release/`, `release-omp/`, `bin/`, or `jobs/`.
- Be careful with `src/wp_def_*.f90`: the number of particles is compiled in, and `build.bash` may edit these files temporarily during builds.
- Preserve scientific behavior unless the user explicitly asks for a change. Numerical code changes should be validated with representative sample inputs when practical.
- Treat `matelem.f90` and `linalg.f90` as performance-sensitive code. Preserve the current algorithmic complexity, no-gradient fast paths, cache-aware loops, and MPI routing unless a task explicitly calls for changing them. For performance changes, check numerical agreement up to expected roundoff and benchmark a representative basis size and process count when practical.
- Matrix-element routine names now identify both their purpose and ECG basis (for example, `MatrixElementsHS_RG_0S` and `MatrixElementsAll_RG_0S`). Use the current names from the affected `matelem.f90`; do not restore historical generic names, and update all callers and analogous sibling codes together when renaming an interface.
- Preserve copyright and license banner comments in bundled third-party sources. Keep `THIRD-PARTY-NOTICES.md` consistent if bundled third-party code or its notices change.

## Building

The preferred batch build entry point is the root script:

```bash
./build.bash machine=linux-generic toolchain=systemdefault config=release code=RG_0S nparticles=4 precision=8 linalg=netlib openmp=no
```

Run `./build.bash` from the repository root; invoking it with no arguments prints usage and exits with status 1. It loops over requested toolchains, configurations, codes, particle counts, precisions, linear-algebra choices, and OpenMP settings, then places binaries under:

```text
bin/[<machine>/]<toolchain>/<config>/<CODE>_N<nparticles>_P<precision>_<linalg>[_omp]
```

The `<machine>/` segment is present only for the cluster machines `puma`, `ocelote`, and `elgato`; it is empty for `linux-generic`, `ubuntu-generic`, `irgetas`, `shabyt`, and `muon`. The `<linalg>` suffix is always present. OpenMP builds are written to the same `<config>` directory as the serial ones and are distinguished by a trailing `_omp` in the file name.

`nparticles` is the only required argument. The defaults are `machine=linux-generic`, `toolchain=systemdefault`, `config=debug,release`, `precision=8`, `linalg=netlib`, `openmp=no`, and a `code` list containing all codes except `CG_0S`. Every argument except `machine` accepts a comma-separated list of values. For `precision=10` or `precision=16`, `build.bash` builds only `linalg=netlib` and skips optimized-library choices; `openmp=yes` is skipped for any code other than `RG_0S`, `RG_1P`, `RG_2D`, and `RG_2P`.

Supported machines are `linux-generic` (default), `ubuntu-generic`, `irgetas`, `shabyt`, `muon`, `puma`, `ocelote`, and `elgato`. Toolchains are machine-specific: the generic and NU machines accept `systemdefault`, `foss-2023b`/`foss-2024a`/`foss-2025a`/`foss-2025b`, `intel-2023b`/`intel-2024a`/`intel-2025a`/`intel-2025b`, and `nvhpc-25.3`/`nvhpc-25.9`/`nvhpc-26.5`; `puma` accepts `systemdefault`, `ohpc`, and `intel-2024.0.0`; `ocelote` and `elgato` accept `systemdefault`, `ohpc`, and `intel-2020.4`. Anything other than `systemdefault` loads an Easybuild or `module` environment. Machine profiles are hardcoded in both `build.bash` and the per-code Makefiles. Toolchain names and module-loading logic live in `build.bash`; the Makefiles select flags by `MACHINE`, `COMPILER`, and precision.

To build a single code directly, use its Makefile:

```bash
cd RG_0S
make release COMPILER=gfortran MACHINE=linux-generic PREC=8 LINALG=openblas EXEFILE=ecg
make release COMPILER=gfortran MACHINE=linux-generic PREC=8 LINALG=openblas OPENMP=yes EXEFILE=ecg
make debug   COMPILER=gfortran MACHINE=linux-generic PREC=8 LINALG=netlib   EXEFILE=ecg
make clean
```

The Makefiles also provide `cleaner`, `cleanest`, `cleanrelease`, and `cleandebug`. Object and module files are written under `release/` or `debug/`; the four real-ECG energy codes use separate `release-omp/` or `debug-omp/` trees when `OPENMP=yes`. Direct Makefile builds place the executable in that same tree using `EXEFILE` unchanged, for example `release/ecg` or `release-omp/ecg`.

For direct Makefile builds, clean before changing compiler, precision, or linear-algebra provider because these choices share an object directory. `build.bash` already runs `make clean` for each build combination.

Common Makefile parameters:

- `CONFIG`: selected by the target, usually `release` or `debug`. `release` optimization flags are compiler-specific (`-O3 -march=native` for gfortran, `-O3 -xhost -fp-model precise` for ifort/ifx, `-O3 -tp=native` for nvfortran); `debug` enables bounds, uninitialized-value, and FPE checking
- `COMPILER`: `gfortran`, `ifort`, `ifx`, or `nvfortran`; these normally map to `mpif90`, `mpiifort`, `mpiifx`, and `mpif90`, respectively
- `MACHINE`: machine profile such as `linux-generic`, `ubuntu-generic`, `irgetas`, `shabyt`, `muon`, `puma`, `ocelote`, or `elgato`
- `PREC`: `8`, `10`, or `16`
- `LINALG`: `netlib`, `mkl`, `lblas`, `openblas`, or `aocl`
- `OPENMP`: `no` (default) or `yes`, supported by the four real-ECG energy codes; enables OpenMP compiler/linker flags and selects a separate build tree
- `EXEFILE`: output executable name for direct Makefile builds

Precision notes:

- `PREC=8`: double precision; broadly supported
- `PREC=10`: extended precision; mainly GNU/x86; use `LINALG=netlib`
- `PREC=16`: quadruple precision; slow emulation on most hardware; use `LINALG=netlib`

The current Makefiles support precisions 8/10/16 with gfortran, 8/16 with ifort or ifx, and only 8 with nvfortran. Compiler support does not guarantee MPI support for every precision; consult `doc/compilation_and_execution.md` before running multi-process high-precision cases.

Linear algebra notes:

- Energy codes can link BLAS/LAPACK.
- `LINALG=netlib` compiles bundled `src/BLAS.f` and `src/LAPACK.f`. These are the single copies used by both the main code and `src/qrupdate/`; the duplicate `qrupdate/BLAS.f` and `qrupdate/LAPACK.f` were removed, so add any new reference routine to the top-level files only.
- Optimized libraries are supported only for `PREC=8`; `PREC=10` and `PREC=16` require `LINALG=netlib`.
- Off-diagonal matrix-element codes accept `LINALG` for interface consistency, but it is a no-op because those codes do not link BLAS/LAPACK.

The number of particles is a compile-time parameter, not a runtime option. `build.bash` accepts particle counts from 2 through 19, temporarily edits `Glob_AllowedNumOfParticles` and `MPI_DPREC` in the selected `src/wp_def_<PREC>.f90`, performs the build, and restores the original file from `wp_def_temporary.f90`. A direct Makefile build uses the values already in that source file. A binary built for one particle count will reject inputs with a different `PARTICLES` value. Avoid concurrent builds of the same code directory, which share the temporary backup and build trees; after an interrupted build, check that the precision-definition file was restored.

## Running

All binaries are MPI programs:

```bash
mpirun -np <NPROCS> /path/to/bin/[<machine>/]<toolchain>/<config>/<binary>
```

Run from a work directory containing the required input files. For validation, copy the sample case into `jobs/` or another scratch directory and use an absolute executable path so calculation output does not overwrite tracked examples.

For OpenMP builds, set `OMP_NUM_THREADS` explicitly and allocate enough cores for the MPI rank count multiplied by the threads per rank. In `LINALG=netlib` builds, OpenMP parallelizes sufficiently large bundled BLAS operations used by qrlinalg. External BLAS providers have their own threading settings; see `doc/compilation_and_execution.md` for placement and threading guidance.

Energy codes read and write a single `inout.txt` in the current working directory; this is the default value of `Glob_DataFileName` in `globvars.f90`. The file begins with an optional `BASIS_TYPE <BASIS>` line and a required `PARTICLES <n>` line.

Off-diagonal matrix-element codes read wave-function files for the two states rather than a single `inout.txt`. Most use `wf_state0.txt` and `wf_state1.txt`; `RG_1P-2D` and `RG_1P-2P` use `wf_state1.txt` and `wf_state2.txt`.

Some fine-structure coupling codes also require a small `inp.txt`:

- no extra `inp.txt`: `RG_0S-1P`, `RG_1P-2D`, `RG_1P-2P`, `RG_0S-2D`
- requires `inp.txt`: `RG_0S-2P`, `RG_1P-1P`, `RG_2D-2D`, `RG_2P-2D`, `RG_2P-2P`

Where present, the positional value or values in `inp.txt` select the requested transition or coupling type.

There is no general automated test suite. Validation is usually done by running sample physical cases and comparing energies, expectation values, or matrix elements against documented reference values.

Sample inputs live under each code's `sample_input/` directory when available. The four real-ECG energy codes and all nine off-diagonal codes provide worked examples; `CG_0S` currently has none. The example READMEs document the calculation, expected runtime, and reference values:

- Energy-code cases are single directories containing `inout.txt`. Recurring case types include `basis_generation_*`, `expected_values_*`, `densities_*`, and `store_wavefunction_*`.
- Transition-dipole cases are named `transition_dipole_moment_*` and put a `README.md` in each case directory.
- Fine-structure cases put instructions in their `initial_state_*/`, `final_state_*/`, and `matelem/` subdirectories rather than at case level. The `matelem/` directory holds both wave-function files and, except for `RG_0S-2D`, the coupling-selector `inp.txt`.
- `RG_0S-1P` and `RG_1P-2D` additionally have an overview `sample_input/README.md`.

The density examples cover both coordinate space (`DENSITIES`) and momentum space (`MOMT_DENS`). They include Python grid builders (`create_grid.py` for the radial `RG_0S` example; `create_cyl_grid.py` and the alternative `create_radial_grid.py` for `RG_1P`, `RG_2D`, and `RG_2P`), supplied grid `.dat` files, and gnuplot scripts `create_plots_{coor,mom}_{dens,cf}.gp`. Reference tables and plots are committed under `output_cyl/` and `cyl_plots/`, including the radial `RG_0S` case despite those directory names. Follow each example's README for grid generation and plotting.

## Source Architecture

In the four real-ECG energy codes, which have the fullest source set, module order is:

```text
wp_def_<PREC> -> globvars -> misc, linalg, spin -> matelem -> matform -> qrupdate -> qrlinalg -> workproc -> main
```

Other codes carry only a subset of these files:

- `CG_0S`: the same set minus `spin.f90`, `qrlinalg.f90`, and `qrupdate/`. It also has no `Q` eigenvalue method and no `MOMT_DENS` or `SAVE_HSWF` step, so do not assume feature parity with the `RG_*` energy codes.
- Off-diagonal codes: only `wp_def_<PREC>.f90`, `globvars.f90`, `matelem.f90`, `workproc.f90`, `main.f90`, and `spin.f90` where spin coupling is needed. They have no `misc`, `linalg`, `matform`, or `qrupdate`, and link no BLAS/LAPACK. `spin.f90` is present in `RG_0S-2D`, `RG_0S-2P`, `RG_1P-1P`, `RG_2D-2D`, `RG_2P-2D`, and `RG_2P-2P`, and absent in the transition-dipole codes `RG_0S-1P`, `RG_1P-2D`, and `RG_1P-2P`.

Important files:

- `wp_def_<PREC>.f90`: working precision kind, MPI real type, and `Glob_AllowedNumOfParticles`
- `globvars.f90`: global state, physical constants, numerical constants, basis/matrix storage; this includes both `Glob_PseudoChargeMatrix` and the pre-scaled `Glob_ScaledPseudoChargeMatrix`
- `misc.f90`: miscellaneous helper routines
- `linalg.f90`: linear-algebra wrappers over BLAS/LAPACK, including performance-aware serial/MPI routing and cache-blocked LDL factorization paths
- `BLAS.f`, `LAPACK.f`: bundled modified reference sources used with `LINALG=netlib`
- `dmng.f`, `X1MACH.f90`: nonlinear optimizer support and machine constants
- `spin.f90`: spin algebra and permutation-symmetry projection
- `matelem.f90`: matrix elements between individual basis functions
- `matform.f90`: assembly of Hamiltonian and overlap matrices
- `qrlinalg.f90` and `qrupdate/`: QR factorization/update state used by method `Q` in the four real-ECG energy codes
- `workproc.f90`: the bulk of the program (roughly 13,000-14,000 lines in the real-ECG energy codes, about 9,000 in `CG_0S`, and 1,300-2,200 in the off-diagonal codes), including `ReadIOFile`/`SaveResults` I/O, basis construction, optimization cycles, generalized symmetric eigensolvers (methods `G`, `I`, and `Q` in the real-ECG energy codes), expectation values, densities, and swap-file handling
- `main.f90`: MPI initialization, random-number seeding, and top-level execution of BBOP input steps

Common BBOP steps handled from `main.f90` include:

- `BASIS_ENL`
- `OPT_CYCLE`
- `FULL_OPT1`
- `ELIM_LCFN`
- `ELIM_LND1`
- `SEPR_LND1`
- `SEPR_FLCF`
- `EXPC_VALS`
- `DENSITIES`
- `MOMT_DENS`
- `SAVE_FILE`
- `SAVE_HSWF`

All four real-ECG energy codes implement this full list. `CG_0S` implements the same steps except `MOMT_DENS` and `SAVE_HSWF`. Adding a new calculation mode usually requires adding a case in `main.f90` and corresponding implementation in `workproc.f90`.

### Current Matrix-Element And Linear-Algebra Conventions

The primary energy-code entry points in `matelem.f90` use symmetry-qualified names:

- Hamiltonian/overlap, optionally with nonlinear-parameter derivatives: `MatrixElementsHS_RG_0S`, `MatrixElementsHS_RG_1P`, `MatrixElementsHS_RG_2D`, `MatrixElementsHS_RG_2P`, and `MatrixElementsHS_CG_0S`
- full expectation-value/operator sets: the corresponding `MatrixElementsAll_*` routines
- normalized or diagonal overlaps, where present: `OverlapMatrixElement_<BASIS>` or `NormalizedOverlapMatElem_<BASIS>`

Off-diagonal codes use specialized operator routines and basis-specific overlap routines. Consult the local `matelem.f90` rather than relying on an older name. Several matrix-element kernels now replace deeply nested particle-pair accumulations with preassembled matrix expressions and skip derivative work when neither gradient is requested. Do not expand these back into the older pair-by-pair algorithms.

`ProgramDataInit` constructs `Glob_ScaledPseudoChargeMatrix` once; hot matrix-element loops should read it directly rather than recomputing scaled charge products. In the energy codes, `linalg_setparam(n)` selects and, for multi-process double-precision runs, calibrates the current linear-algebra paths. It must remain a collective call on all MPI ranks and must be called again when the problem dimension changes.

For method `Q`, preserve canonical basis order and the authoritative lower triangles of `Glob_H` and `Glob_S`, which hold normalized, unshifted physical matrix elements. Upper triangles are scratch or unspecified storage. Consult `doc/devnotes/Q_method_design.md` before changing Q matrix assembly, updates, or solver state.

## Documentation To Check

- Build and execution details: `doc/compilation_and_execution.md`
- Input format: `doc/input_file_format.md`
- Physics and notation: `doc/theoretical_background.md`
- VS Code setup: `doc/use_of_visual_studio_code.md`
- Q-method architecture, canonical H/S layout, porting steps, and validation: `doc/devnotes/Q_method_design.md`

For behavior questions, prefer these project docs over inference from code alone.

## Editing And Validation Checklist

Before editing:

- Identify which code directories are affected.
- Check whether the same pattern appears in sibling `RG_*` directories.
- Review the relevant project documentation if touching input files, build behavior, or physics conventions.

After editing:

- Run a focused build when possible, usually with `make debug` or `make release` for the affected code.
- For numerical behavior changes, run a small sample input if practical.
- Check for unintended edits to generated files, temporary files, `debug/`, `release/`, `bin/`, and `jobs/`.
- Keep `AGENTS.md` and `CLAUDE.md` conceptually aligned when changing agent-facing project guidance.
