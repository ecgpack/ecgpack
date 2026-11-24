# AGENTS.md

Guidance for AI coding agents working in this repository. This file is the single source of agent instructions: Codex and Antigravity read it directly, and `CLAUDE.md` imports it for Claude Code.

## Project Overview

ECGPACK is a collection of closely related parallel Fortran codes for high-accuracy variational calculations of quantum few-body systems, including few-electron atoms, molecules, ions, and systems with exotic particles. The codes use all-particle explicitly correlated Gaussian (ECG) basis functions and MPI parallelism. The public repository is <https://github.com/ecgpack/ecgpack>, and its default branch is `main`.

The project is developed by the research group of Sergiy Bubin, Physics Department, Nazarbayev University. For scientific background, notation, and references, see `README.md` and the documents in `doc/`.

## Repository Layout

Each main code lives in a top-level directory with a similar structure:

- `Makefile`
- `src/`
- `.vscode/`: per-code VS Code tasks and launch configurations; see `doc/use_of_visual_studio_code.md`
- usually `sample_input/`

The main code groups overlap:

- Energy and wave-function codes: `RG_0S`, `RG_1P`, `RG_2D`, `RG_2P`, `PG_0S`, `PG_1P`, `CG_0S`. These generate ECG basis sets for states of different angular momentum and parity, and can compute energies, expectation values, particle distributions, and saved wave functions.
- Off-diagonal matrix-element codes: `RG_0S-1P`, `RG_0S-2D`, `RG_0S-2P`, `RG_1P-1P`, `RG_1P-2D`, `RG_1P-2P`, `RG_2D-2D`, `RG_2P-2D`, `RG_2P-2P`. These work with pairs of states expanded in separate symmetry-adapted bases; the states may have the same angular momentum and parity. `RG_1P-1P`, `RG_2D-2D`, and `RG_2P-2P` support scalar as well as spin-dependent operators.
- Transition dipole codes: `RG_0S-1P`, `RG_1P-2D`, `RG_1P-2P`. These evaluate transition electric-dipole moments in the length and velocity gauges.
- Fine-structure coupling codes: `RG_0S-2D`, `RG_0S-2P`, `RG_1P-1P`, `RG_2D-2D`, `RG_2P-2D`, `RG_2P-2P`. These evaluate spin-orbit and non-contact spin-spin matrix elements.
- Real ECG codes: all `RG_*` directories.
- Power Gaussian codes: `PG_0S` ($L=0$, even parity) and `PG_1P` ($L=1$, odd parity). These use real ECGs premultiplied by an even power (up to 250) of $r_1$, the distance between particles 1 and 2. They are intended for systems where particles 1 and 2 are heavy and carry charges of the same sign, such as diatomic molecules.
- Complex ECG codes: `CG_0S` (complex $L=0$ basis). It mirrors `RG_0S` (including its CUDA backend), but has no semi-relativistic kinetic energy.

`RG_0S` is the most complete reference implementation. Other codes intentionally share much of its source layout and Makefile structure.

Non-code directories:

- `doc/`: project documentation (`compilation_and_execution.md`, `input_file_format.md`, `theoretical_background.md`, `use_of_visual_studio_code.md`); `doc/devnotes/` contains developer-facing design and maintenance guides (`Q_method_design.md`, the `qrlinalg.f90` generator and test in `qrlinalg_generator/`, and the high-precision Bessel-function notes in `mod_Bessel_high_precision/`), and `doc/img/` holds the project logos
- `utilities/`: helper scripts — input-file generation and merging (`input_file_manipulation/`), nuclear masses (`nuclear_masses/`), job monitoring (`job_monitoring/`), and a Young-tableaux spreadsheet (`permutational_symmetry/`)
- `.agents/skills/`: skills shared by all AI coding agents, with identical copies for Claude Code in `.claude/skills/`; see Agent Skills below
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
- `CLAUDE.md`: a stub that imports this file for Claude Code (`@AGENTS.md`); do not add guidance there

## General Agent Rules

- Prefer small, targeted changes that follow the existing code style.
- Do not rewrite broad portions of the duplicated Fortran code unless the task explicitly requires it.
- Treat the `RG_*` and `PG_*` directories as related implementations. When changing shared algorithms, matrix elements, build logic, or file formats, check whether analogous changes are needed in sibling directories.
- Do not commit generated build products from `debug/`, `debug-omp/`, `release/`, `release-omp/`, `bin/`, or `jobs/`.
- Text files use LF line endings; do not introduce CRLF.
- Be careful with `src/wp_def_*.f90`: the number of particles is compiled in, and `build.bash` edits these files temporarily during builds.
- Preserve scientific behavior unless the user explicitly asks for a change. Numerical code changes should be validated with representative sample inputs when practical.
- Treat `matelem.f90`, `matform.f90`, `linalg.f90`, and their callers in `workproc.f90` as performance-sensitive code. Preserve the current algorithmic complexity, precomputation outside pair/permutation loops, no-gradient fast paths, cache-aware loops, and MPI routing unless a task explicitly calls for changing them. For performance changes, check numerical agreement up to expected roundoff and benchmark a representative basis size and process count when practical.
- Matrix-element routine names identify both their purpose and ECG basis (for example, `MatrixElementsHS_RG_0S`; the PG codes use `MatrixElementsOpt` or `MatrixElements`). Use the current names from the affected `matelem.f90`; do not restore historical generic names, and update all callers and analogous sibling codes together when renaming an interface.
- Preserve copyright and license banner comments in bundled third-party sources. Keep `THIRD-PARTY-NOTICES.md` consistent if bundled third-party code or its notices change.
- Keep this file small. Codex truncates project instructions beyond 32 KiB and Antigravity truncates each instruction file beyond 24,000 bytes, so keep `AGENTS.md` below 22 KB: put detailed reference material in `doc/devnotes/` and link to it from here.
- Agent-specific rule directories (`.claude/rules/`, `.agents/rules/`) are intentionally not used (their formats are not portable); project rules belong in this file.

## Building

Full details (machines, toolchains, precision and linear-algebra support, OpenMP, CUDA, runtime switches) are in `doc/compilation_and_execution.md`; `./build.bash` with no arguments prints its usage and exits with status 1. The essentials:

```bash
./build.bash machine=linux-generic toolchain=systemdefault config=release code=RG_0S nparticles=4 precision=8 linalg=netlib openmp=no
```

- Run it from the repository root. Only `nparticles` (2 through 19) is required. Defaults: `machine=linux-generic`, `toolchain=systemdefault`, `config=debug,release`, `precision=8`, `linalg=netlib`, `openmp=no`, `gpu=no`, and every code `build.bash` knows. Every argument except `machine`, `gpu`, and `cuda_arch` accepts a comma-separated list.
- Binaries go to `bin/[<machine>/]<toolchain>/<config>/<CODE>_N<nparticles>_P<precision>_<linalg>[_omp][_gpu]`. The `<machine>/` segment appears only for `puma`, `ocelote`, and `elgato`; the `<linalg>` suffix is always present.
- `build.bash` does not know `PG_1P` yet (it rejects `code=PG_1P`); build that code with its Makefile.
- Skipped or rejected combinations: `precision=10`/`16` build only `linalg=netlib`; `openmp=yes` applies only to `RG_0S`, `RG_1P`, `RG_2D`, and `RG_2P`; `gpu=yes` requires an `nvhpc-*` toolchain, `precision=8`, and one of `CG_0S`, `PG_0S`, `RG_0S`, `RG_1P`, `RG_2D`, `RG_2P`.
- Machine profiles are hardcoded in both `build.bash` and the per-code Makefiles (adding one means editing both); toolchain names and module loading live only in `build.bash`.

To build a single code directly, use its Makefile:

```bash
cd RG_0S
make release COMPILER=gfortran MACHINE=linux-generic PREC=8 LINALG=openblas EXEFILE=ecg
make debug   COMPILER=gfortran MACHINE=linux-generic PREC=8 LINALG=netlib   EXEFILE=ecg
make clean   # also: cleaner, cleanest, cleanrelease, cleandebug
```

- Parameters: `COMPILER` (`gfortran`, `ifort`, `ifx`, or `nvfortran`; the matching MPI wrapper is what is actually called), `MACHINE`, `PREC` (`8`, `10`, `16`), `LINALG` (`netlib`, `mkl`, `lblas`, `openblas`, `aocl`), `OPENMP` (`no`/`yes`), `USE_CUDA` (`no`/`yes`, with `CUDA_ARCH=sm_XX`), and `EXEFILE`.
- Objects, `.mod` files, and the executable (named `EXEFILE` unchanged) go to `release/` or `debug/`, or to `release-omp/` or `debug-omp/` with `OPENMP=yes`.
- `COMPILER`, `PREC`, `LINALG`, and `USE_CUDA` share one object directory: run `make clean` before changing any of them. `build.bash` already does this for each build combination.
- Parallel `make -j` is safe because the Makefile rules list every module dependency (for example, `spin.o` depends on `globvars.o`). When adding a source file or a `use` statement, add the matching prerequisite to the Makefile rule; otherwise `-j` builds can fail with a missing `.mod` file.
- Supported precisions are 8/10/16 with gfortran, 8/16 with ifort or ifx, and 8 with nvfortran. Compile support is not run support: gfortran with `PREC=16` currently runs with a single MPI process only.
- `LINALG=netlib` compiles the bundled `src/BLAS.f` and `src/LAPACK.f`, the single copies used by both the main code and `src/qrupdate/`; add any new reference routine to these files only. Optimized libraries require `PREC=8`. In the off-diagonal codes `LINALG` is accepted but is a no-op because they link no BLAS/LAPACK.
- `OPENMP=yes` is supported by the four real-ECG energy codes only (the `PG_*` Makefiles accept it, but OpenMP builds of the PG codes are not supported). The `!$OMP` directives live in those codes' `src/BLAS.f` and are gated by an `OMP_MIN_WORK` size threshold.
- `USE_CUDA=yes` compiles the CUDA Fortran backend `src/gpu_backend.f90` in `RG_0S`, `RG_1P`, `RG_2D`, `RG_2P`, `PG_0S`, and `CG_0S`; it requires `COMPILER=nvfortran` and `PREC=8`. The binary still computes on the CPU unless `ECG_GPU=1` or `ECG_GPU_EIG=1` is set at runtime, and it needs a visible CUDA device either way.

The number of particles is a compile-time parameter, not a runtime option. `build.bash` temporarily edits `Glob_AllowedNumOfParticles` and `MPI_DPREC` in the selected `src/wp_def_<PREC>.f90`, performs the build, and restores the original file from `wp_def_temporary.f90`; a direct Makefile build uses the values already in that source file. A binary built for one particle count rejects inputs (or, for the off-diagonal codes, wave-function files) with a different `PARTICLES` value; every code requires an exact match. Avoid concurrent builds of the same code directory (they share the backup file and build trees), and after an interrupted build check that `src/wp_def_<PREC>.f90` was restored.

## Running

All binaries are MPI programs:

```bash
mpirun -np <NPROCS> /path/to/bin/[<machine>/]<toolchain>/<config>/<binary>
```

Run from a work directory containing the required input files. For validation, copy the sample case into `jobs/` or another scratch directory and use an absolute executable path so calculation output does not overwrite tracked examples.

For OpenMP builds, set `OMP_NUM_THREADS` explicitly and allocate enough cores for the MPI rank count multiplied by the threads per rank. External BLAS providers have their own threading settings; see `doc/compilation_and_execution.md` for placement and threading guidance.

Energy codes read and write a single `inout.txt` in the current working directory; this is the default value of `Glob_DataFileName` in `globvars.f90`. The file begins with an optional `BASIS_TYPE <BASIS>` line and a required `PARTICLES <n>` line. Some header keywords are basis-specific: `FIXED_INDEX` (`RG_1P`, `PG_1P`), `FIXED_INDEX_1`/`FIXED_INDEX_2` (`RG_2D`, `RG_2P`), and `VECTOR_COUPLING_SCHEME` (`RG_2D` only); see `doc/input_file_format.md`.

Off-diagonal matrix-element codes read wave-function files for the two states rather than a single `inout.txt`. Most use `wf_state0.txt` and `wf_state1.txt`; `RG_1P-2D` and `RG_1P-2P` use `wf_state1.txt` and `wf_state2.txt`.

Some fine-structure coupling codes also require a small `inp.txt`:

- no extra `inp.txt`: `RG_0S-1P`, `RG_1P-2D`, `RG_1P-2P`, `RG_0S-2D`
- requires `inp.txt`: `RG_0S-2P`, `RG_1P-1P`, `RG_2D-2D`, `RG_2P-2D`, `RG_2P-2P`

Where present, the positional value or values in `inp.txt` select the requested transition or coupling type. In `RG_1P-1P`, `RG_2D-2D`, and `RG_2P-2P`, the first line contains two tokens: a mode (`SPIN`, `SCALAR`, or `ALL`) and a transition (`XP_XP`/`3P_1P` for the P codes, `XD_XD`/`3D_1D` for the D code). The triplet-to-singlet transitions accept only `SPIN`. Consult each code's `main.f90` and sample README for its selectors; they are not interchangeable across codes.

There is no general automated test suite. Validation is usually done by running sample physical cases and comparing energies, expectation values, or matrix elements against documented reference values.

Sample inputs live under each code's `sample_input/` directory; every code provides worked examples. The example READMEs document the calculation, expected runtime, and reference values:

- Energy-code cases are single directories containing `inout.txt`. Recurring case types include `basis_generation_*`, `expected_values_*`, `densities_*`, and `store_wavefunction_*`. `PG_0S` adds `basis_health_checking_*` (`SAVE_HS_R`), `linear_c_elimination_*` (`ELIM_LCFN`), `linear_depen_elimination_*` (`ELIM_LND1`), and `overlap_diagonalization_*` (`OVERLAP_D`).
- Transition-dipole cases are named `transition_dipole_moment_*` and put a `README.md` in each case directory.
- Fine-structure cases put instructions in their `initial_state_*/`, `final_state_*/`, and `matelem/` subdirectories rather than at case level. The `matelem/` directory holds both wave-function files and, except for `RG_0S-2D`, the coupling-selector `inp.txt`.
- `RG_0S-1P` and `RG_1P-2D` additionally have an overview `sample_input/README.md`.

The `RG_*` density examples (`DENSITIES`, `MOMT_DENS`) ship grid builders, gnuplot scripts, and reference output; the PG ones compute only the nucleus–nucleus correlation function. Follow each example's README.

The `run-sample-case` skill automates this for the energy codes (see Agent Skills below).

## Utilities

- `utilities/input_file_manipulation/ecgpack_input_file_generator.py`: generates an initial input file on standard output; use `python3 ... -h` for its options and redirect output into a job directory. Its `-particles` option currently accepts 2 through 9, a narrower range than `build.bash`, and its `-basis_type` option accepts the four `RG_*` energy bases, `PG_0S`, `PG_1P`, and `CG_0S`.
- `utilities/input_file_manipulation/ecgpack_input_file_merger.py`: concatenates selected basis ranges from two inputs using `-file1`, `-file2`, and `-file3`, with optional `-range1` and `-range2`. It takes header settings from the first file and writes an `EXPC_VALS G` command for the merged basis. Preserve its verbatim copying of numerical strings so extended and quadruple precision are not lost through Python float conversion. It does not check that the two inputs describe the same system (masses, charges, symmetry).
- `utilities/nuclear_masses/nuclear_mass_determination.py`: prints nuclear masses and uncertainty estimates using embedded CODATA 2022 constants and AME 2020 atomic masses; retain the source references when updating its data.
- `utilities/job_monitoring/status_all.bash`: reports calculation progress and optional SLURM status. Its `JOBSDIR`, directory selection, and display flags are configured near the top; review these local settings before running it.
- `utilities/permutational_symmetry/Young_Tableaux.xlsx`: spreadsheet supporting permutation-symmetry work.

## Source Architecture

There is no separate architecture document: this section and the comments in the sources are the reference. **Read the header comments of the routines you touch before editing `globvars.f90`, `matelem.f90`, `matform.f90`, `linalg.f90`, `qrlinalg.f90`, `workproc.f90`, `main.f90`, or `gpu_backend.f90`.** In brief:

- Module order in the seven energy codes: `wp_def_<PREC> -> globvars -> misc, linalg, spin -> matelem -> matform -> qrlinalg -> workproc -> main` (`qrlinalg` needs only `wp_def`; `src/qrupdate/` is compiled and linked but unused). CUDA builds add `gpu_backend` after `matelem`. The off-diagonal codes carry only `wp_def_<PREC>`, `globvars`, `matelem`, `workproc`, `main`, and `spin` where spin coupling is needed.
- `main.f90` executes a sequence of BBOP (Basis Building and Optimization Program) steps read from the input: `BASIS_ENL`, `OPT_CYCLE`, `FULL_OPT1`, `ELIM_LCFN`, `ELIM_LND1`, `SEPR_LND1`, `SEPR_FLCF`, `EXPC_VALS`, `DENSITIES`, `MOMT_DENS`, `SAVE_FILE`, `SAVE_HSWF`. The PG codes lack `MOMT_DENS` and add `SAVE_HS_R` and `OVERLAP_D`. Adding a calculation mode means adding a `case` in `main.f90` and a routine in `workproc.f90`, which holds the bulk of each program (I/O, basis enlargement, optimization, the eigensolvers `G`, `I`, and `Q`, expectation values, densities).

Invariants that are easy to break:

- `Glob_n`, `Glob_np`, `Glob_npt`, `Glob_2Raised3n2`, and `Glob_PiRaised3n2` are named constants fixed by the compiled-in particle count. Never assign or `MPI_BCAST` them, and pass full-size arrays to explicit-shape dummy arguments even on paths that never touch them (for example `Dk(2*Glob_np)` for `MatrixElementsHS_*`).
- Hot paths read the precomputed `Glob_ScaledPseudoChargeMatrix` and the once-per-sweep outputs of `PrecomputeMatrices_L_A` (`PrecomputeMatrices_MA_L_A` in `PG_0S`, `PrecomputeMatrices_L_C` in `CG_0S`); do not recompute them per basis-function pair.
- `linalg_setparam(n)` must remain a collective call on all MPI ranks and be repeated whenever the problem dimension changes.
- Method `Q`: basis arrays stay in canonical order, and the unshifted lower triangles of `Glob_H` and `Glob_S` are authoritative. Read `doc/devnotes/Q_method_design.md` first. `qrlinalg.f90` is identical in all seven energy codes, and every call into it is collective: never guard one by rank or broadcast its results. Its complex half (for `CG_0S`) is generated; never edit it by hand. The Q helpers in `workproc.f90` are the `RG_0S` text in every code (PG style; Hermitian in `CG_0S`), so a fix goes into all seven.
- `misc.f90` and the RNG-seeding `block` in `main.f90` are byte-identical in all seven energy codes; keep them that way. Use `urand` and `drnor` from `misc.f90`, never the intrinsic `random_number`/`random_seed`.
- Each MPI rank seeds its generators differently, so a random choice that all ranks must agree on has to be made on rank 0 and broadcast, as `GenerateRndIntSeq` does.
- CUDA builds share one arithmetic kernel between CPU and GPU (`MatrixElementsHS_<BASIS>`, or `MatrixElementsOpt` in `PG_0S`, marked `attributes(host,device)`). Do not duplicate the arithmetic, and keep the NVHPC workarounds noted in the comments of those kernels and of `gpu_backend.f90` (for example, no I/O in device code, and `cudaMemcpy` instead of array-section device assignments).

## Agent Skills

Skills shared by all agents live in `.agents/skills/<skill-name>/` (a `SKILL.md` plus optional `scripts/`), the location that Codex and Antigravity read natively. Claude Code reads `.claude/skills/` instead, so each skill is also kept there as an identical copy (not a symlink, which can break in Windows checkouts). Edit a skill in `.agents/skills/`, then refresh the copy:

```bash
rm -rf .claude/skills/<skill-name> && cp -r .agents/skills/<skill-name> .claude/skills/
diff -r .agents/skills .claude/skills   # must print nothing
```

Refer to helper scripts by their `.agents/skills/...` path so that both copies stay identical. Restrict the `SKILL.md` frontmatter to the fields of the Agent Skills standard (`name`, which must match the directory name, and `description`; optionally `license`, `compatibility`, and `metadata`) so that every agent accepts it.

Available skills:

- `run-sample-case`: builds an energy code if needed and runs one of its `sample_input` cases in a scratch copy under `jobs/`.

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
- Check whether the same pattern appears in sibling `RG_*` or `PG_*` directories.
- Review the relevant project documentation if touching input files, build behavior, or physics conventions.

After editing:

- Run a focused build when possible, usually with `make debug` or `make release` for the affected code.
- When changing module declarations or interfaces, also build with an Intel toolchain if one is available: gfortran accepts some non-standard constructs that ifx rejects, such as an entity given the `PUBLIC` attribute twice.
- For numerical behavior changes, run a small sample input if practical.
- Check for unintended edits to generated files, temporary files, `debug/`, `release/`, `bin/`, and `jobs/`.
