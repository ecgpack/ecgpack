# AGENTS.md

Guidance for AI coding agents working in this repository.

## Project Overview

ECGPACK is a collection of closely related parallel Fortran codes for high-accuracy variational calculations of quantum few-body systems, including few-electron atoms, molecules, ions, and systems with exotic particles. The codes use all-particle explicitly correlated Gaussian (ECG) basis functions and MPI parallelism.

The project is developed by the research group of Sergiy Bubin, Physics Department, Nazarbayev University. For scientific background, notation, and references, see `README.md` and the documents in `doc/`.

## Repository Layout

Each main code lives in a top-level directory with a similar structure:

- `Makefile`
- `src/`
- `.vscode/`
- usually `sample_input/`

The main code groups overlap:

- Energy and wave-function codes: `RG_0S`, `RG_1P`, `RG_2D`, `RG_2P`, `CG_0S`. These generate ECG basis sets for states of different angular momentum and parity, and can compute energies, expectation values, particle distributions, and saved wave functions.
- Off-diagonal matrix-element codes: `RG_0S-1P`, `RG_0S-2D`, `RG_0S-2P`, `RG_1P-1P`, `RG_1P-2D`, `RG_1P-2P`, `RG_2D-2D`, `RG_2P-2D`, `RG_2P-2P`. These work with pairs of states expanded in different symmetry-adapted bases.
- Transition dipole codes: `RG_0S-1P`, `RG_1P-2D`, `RG_1P-2P`. These evaluate transition electric-dipole moments in the length and velocity gauges.
- Fine-structure coupling codes: `RG_0S-2D`, `RG_0S-2P`, `RG_1P-1P`, `RG_2D-2D`, `RG_2P-2D`, `RG_2P-2P`. These evaluate spin-orbit and non-contact spin-spin matrix elements.
- Real ECG codes: all `RG_*` directories.
- Complex ECG codes: `CG_0S`. This code currently supports only the complex $L=0$ basis and is a work in progress that may lag behind the real-ECG codes.

`RG_0S` is the most complete reference implementation. Other codes intentionally share much of its source layout and Makefile structure.

Non-code directories:

- `doc/`: project documentation
- `utilities/`: utility scripts and related tools
- `bin/`: generated/user-created binaries, not normally committed
- `jobs/`: user-created calculation work directories, not normally committed

Root files of interest:

- `README.md`: main project manual
- `AUTHORS.md`: contributor list
- `CITATION.cff`: machine-readable citation metadata
- `LICENSE.md`: BSD 3-Clause license for ECGPACK
- `THIRD-PARTY-NOTICES.md`: notices for bundled BLAS, LAPACK, PORT, and SLATEC-derived sources
- `build.bash`: batch build driver
- `.code-workspace`: VS Code multi-folder workspace
- `CLAUDE.md`: Claude-specific project context; keep it consistent with this file when relevant

## General Agent Rules

- Prefer small, targeted changes that follow the existing code style.
- Do not rewrite broad portions of the duplicated Fortran code unless the task explicitly requires it.
- Treat the `RG_*` directories as related implementations. When changing shared algorithms, matrix elements, build logic, or file formats, check whether analogous changes are needed in sibling directories.
- Do not commit generated build products from `debug/`, `release/`, `bin/`, or `jobs/`.
- Be careful with `src/wp_def_*.f90`: the number of particles is compiled in, and `build.bash` may edit these files temporarily during builds.
- Preserve scientific behavior unless the user explicitly asks for a change. Numerical code changes should be validated with representative sample inputs when practical.
- Treat `matelem.f90` and `linalg.f90` as performance-sensitive code. Preserve the current algorithmic complexity, no-gradient fast paths, cache-aware loops, and MPI routing unless a task explicitly calls for changing them. For performance changes, check numerical agreement up to expected roundoff and benchmark a representative basis size and process count when practical.
- Matrix-element routine names now identify both their purpose and ECG basis (for example, `MatrixElementsHS_RG_0S` and `MatrixElementsAll_RG_0S`). Use the current names from the affected `matelem.f90`; do not restore historical generic names, and update all callers and analogous sibling codes together when renaming an interface.
- Preserve copyright and license banner comments in bundled third-party sources. Keep `THIRD-PARTY-NOTICES.md` consistent if bundled third-party code or its notices change.

## Building

The preferred batch build entry point is the root script:

```bash
./build.bash machine=linux-generic toolchain=systemdefault config=release code=RG_0S nparticles=4 precision=8 linalg=netlib
```

Run `./build.bash` with no arguments to see its usage. It loops over requested toolchains, configurations, codes, particle counts, precisions, and linear-algebra choices, then places binaries under:

```text
bin/<toolchain>/<config>/<CODE>_N<nparticles>_P<precision>_<linalg>
```

`nparticles` is the only required argument. The other arguments have defaults; notably, the default precision is `8` and the default linear-algebra choice is `netlib`. For `precision=10` or `precision=16`, `build.bash` builds only `linalg=netlib` and skips optimized-library choices.

To build a single code directly, use its Makefile:

```bash
cd RG_0S
make release COMPILER=gfortran MACHINE=linux-generic PREC=8 LINALG=openblas EXEFILE=ecg
make debug   COMPILER=gfortran MACHINE=linux-generic PREC=8 LINALG=netlib   EXEFILE=ecg
make clean
```

The Makefiles also provide `cleaner`, `cleanest`, `cleanrelease`, and `cleandebug`. Object and module files are written under the selected `release/` or `debug/` directory.

Common Makefile parameters:

- `CONFIG`: selected by the target, usually `release` or `debug`
- `COMPILER`: `gfortran`, `ifort`, `ifx`, or `nvfortran`; these normally map to `mpif90`, `mpiifort`, `mpiifx`, and `mpif90`, respectively
- `MACHINE`: machine profile such as `linux-generic`, `ubuntu-generic`, `irgetas`, `shabyt`, `muon`, `puma`, `ocelote`, or `elgato`
- `PREC`: `8`, `10`, or `16`
- `LINALG`: `netlib`, `mkl`, `lblas`, `openblas`, or `aocl`
- `EXEFILE`: output executable name for direct Makefile builds

Precision notes:

- `PREC=8`: double precision; broadly supported
- `PREC=10`: extended precision; mainly GNU/x86; use `LINALG=netlib`
- `PREC=16`: quadruple precision; slow emulation on most hardware; use `LINALG=netlib`

Linear algebra notes:

- Energy codes can link BLAS/LAPACK.
- `LINALG=netlib` compiles bundled `src/BLAS.f` and `src/LAPACK.f`.
- Optimized libraries are supported only for `PREC=8`; `PREC=10` and `PREC=16` require `LINALG=netlib`.
- Off-diagonal matrix-element codes accept `LINALG` for interface consistency, but it is a no-op because those codes do not link BLAS/LAPACK.

The number of particles is a compile-time parameter, not a runtime option. `build.bash` temporarily edits `Glob_AllowedNumOfParticles` in the selected `src/wp_def_<PREC>.f90`, performs the build, and restores the original file from `wp_def_temporary.f90`. A binary built for one particle count will reject inputs with a different `PARTICLES` value.

## Running

All binaries are MPI programs:

```bash
mpirun -np <NPROCS> /path/to/bin/<toolchain>/<config>/<binary>
```

Run from a work directory containing the required input files.

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

## Source Architecture

Within each `src/`, module order is generally:

```text
wp_def_<PREC> -> globvars -> misc, linalg, spin -> matelem -> matform -> workproc -> main
```

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
- `workproc.f90`: the bulk of the program, including `ReadIOFile`/`SaveResults` I/O, basis construction, optimization cycles, generalized symmetric eigensolvers (methods `G` and `I`), expectation values, densities, and swap-file handling
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

Adding a new calculation mode usually requires adding a case in `main.f90` and corresponding implementation in `workproc.f90`.

### Current Matrix-Element And Linear-Algebra Conventions

The primary energy-code entry points in `matelem.f90` use symmetry-qualified names:

- Hamiltonian/overlap, optionally with nonlinear-parameter derivatives: `MatrixElementsHS_RG_0S`, `MatrixElementsHS_RG_1P`, `MatrixElementsHS_RG_2D`, `MatrixElementsHS_RG_2P`, and `MatrixElementsHS_CG_0S`
- full expectation-value/operator sets: the corresponding `MatrixElementsAll_*` routines
- normalized or diagonal overlaps, where present: `OverlapMatrixElement_<BASIS>` or `NormalizedOverlapMatElem_<BASIS>`

Off-diagonal codes use specialized operator routines and basis-specific overlap routines. Consult the local `matelem.f90` rather than relying on an older name. Several matrix-element kernels now replace deeply nested particle-pair accumulations with preassembled matrix expressions and skip derivative work when neither gradient is requested. Do not expand these back into the older pair-by-pair algorithms.

`ProgramDataInit` constructs `Glob_ScaledPseudoChargeMatrix` once; hot matrix-element loops should read it directly rather than recomputing scaled charge products. In the energy codes, `linalg_setparam(n)` selects and, for multi-process double-precision runs, calibrates the current linear-algebra paths. It must remain a collective call on all MPI ranks and must be called again when the problem dimension changes.

## Documentation To Check

- Build and execution details: `doc/compilation_and_execution.md`
- Input format: `doc/input_file_format.md`
- Physics and notation: `doc/theoretical_background.md`
- VS Code setup: `doc/use_of_visual_studio_code.md`

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
