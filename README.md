<p align="center">
  <img src="doc/img/logo_and_name.svg" alt="ECG package logo" width="480">
</p>

# ECGPACK

ECGPACK is a collection of closely related parallel computer codes for performing high accuracy variational calculations of quantum few-body systems such as few-electron atoms, molecules, ions, and systems containing exotic particles, using all-particle explicitly correlated Gaussian (ECG) basis sets. ECGPACK is written in Fortran and uses MPI for parallelism.

## Retrieving ECGPACK source code [![Latest release](https://img.shields.io/github/v/release/ecgpack/ecgpack)](https://github.com/ecgpack/ecgpack/releases)

To clone the current main (development) branch:

```bash
git clone git@github.com:ecgpack/ecgpack.git ecgpack
```

To clone a release version (e.g. `vXXXX.YY.ZZ`)

```bash
git clone -b vXXXX.YY.ZZ --recursive git@github.com:ecgpack/ecgpack.git ecgpack-vXXXX.YY.ZZ
```

## Citation

Machine-readable citation metadata is provided in [`CITATION.cff`](CITATION.cff) in the root directory.

## Theoretical background

For theoretical background and mathematical notations please see [Theoretical background](doc/theoretical_background.md) in the documetation folder, where relevant references are provided.

## Directory structure

The ECGPACK project repository has the following directory structure:

| Directory | Description |
| --- | --- |
| `ecgpack/` | Root directory |
| `ecgpack/CG_0S/` | The code for energy and wavefunction calculations with complex spherically symmetric Gaussians ($L=0$, even parity) that have the form $$\phi_k=\exp [ \mathbf{r}' (C_k \otimes \mathbf{I}) \mathbf{r}]=\exp [ \mathbf{r}' ((A_k+i B_k) \otimes \mathbf{I}) \mathbf{r}].$$ The CG_0S code is currently a work in progress and lags behind the corresponding RG_0S code in terms of features and the quality of the implementation. |
| `ecgpack/RG_0S/` | The code for energy and wavefunction calculations with real spherically symmetric Gaussians ($L=0$, even parity) that have the form $$\phi_k=\exp [ \mathbf{r}' (A_k \otimes \mathbf{I}) \mathbf{r}].$$ |
| `ecgpack/RG_0S-1P/` | The code for the calculation of the transition electric dipole moments for states that are expanded using RG_0S and RG_1P bases. Currently this includes the transition dipole moment operators in the length and velocity gauges. |
| `ecgpack/RG_0S-2D/` | The code for the calculation of the offdiagonal matrix elements between states that are expanded using RG_0S and RG_2D bases. Currently this includes the evaluation of the spin--orbit and noncontact spin-spin interactions. |
| `ecgpack/RG_0S-2P/` | The code for the calculation of the offdiagonal matrix elements between states that are expanded using RG_0S and RG_2P bases. Currently this includes the evaluation of the spin--orbit and noncontact spin-spin interactions. |
| `ecgpack/RG_1P/` | The code for energy and wavefunction calculations with real ECGs ($L=1$, odd parity) that have the form $$\phi_k=z_{i_k}\exp [ \mathbf{r}' (A_k \otimes \mathbf{I}) \mathbf{r}].$$ |
| `ecgpack/RG_1P-1P/` | The code for the calculation of the offdiagonal matrix elements between states that are expanded using RG_1P and RG_1P bases. Currently this includes the evaluation of the spin--orbit and noncontact spin-spin interactions. |
| `ecgpack/RG_1P-2D/` | The code for the calculation of the transition electric dipole moments for states that are expanded using RG_1P and RG_2D bases. Currently this includes the transition dipole moment operators in the length and velocity gauges. |
| `ecgpack/RG_1P-2P/` | The code for the calculation of the transition electric dipole moments for states that are expanded using RG_1P and RG_2P bases. Currently this includes the transition dipole moment operators in the length and velocity gauges. |
| `ecgpack/RG_2D/` | The code for energy and wavefunction calculations with real ECGs ($L=2$, even parity) that have the form $$\phi_k=(x_{i_k} x_{j_k} + y_{i_k} y_{j_k} - 2 z_{i_k} z_{j_k} ) \exp [ \mathbf{r}' (A_k \otimes \mathbf{I}) \mathbf{r}],$$ where the integer index $i_k$ can be either different or the same as $j_k$. |
| `ecgpack/RG_2D-2D/` | The code for the calculation of the offdiagonal matrix elements between states that are expanded using RG_2D and RG_2D bases. Currently this includes the evaluation of the spin--orbit and noncontact spin-spin interactions. |
| `ecgpack/RG_2P/` | The code for energy and wavefunction calculations with real ECGs ($L=1$, even parity) that have the form $$\phi_k=(x_{i_k} y_{j_k} - x_{j_k} y_{i_k}) \exp [ \mathbf{r}' (A_k \otimes \mathbf{I}) \mathbf{r}].$$ |
| `ecgpack/RG_2P-2D/` | The code for the calculation of the offdiagonal matrix elements between states that are expanded using RG_2P and RG_2D bases. Currently this includes the evaluation of the spin--orbit and noncontact spin-spin interactions. |
| `ecgpack/RG_2P-2P/` | The code for the calculation of the offdiagonal matrix elements between states that are expanded using RG_2P and RG_2P bases. Currently this includes the evaluation of the spin--orbit and noncontact spin-spin interactions. |
| `ecgpack/bin/` | Binary files for calculations (may be created by user) |
| `ecgpack/doc/` | Directory containing manuals and documentation |
| `ecgpack/jobs/` | Work directory for calculations (may be created by user) |
| `ecgpack/utilities/` | Various utilities and scripts |

All codes can be divided into three main groups:

| Group name | List of Codes | Description |
| --- | --- | --- |
| Energy and wavefunctions codes | `CG_0S`, `RG_0S`, `RG_1P`, `RG_2D`, `RG_2P` | Generation of ECG basis sets for states of different angular-momentum/parity. They can compute energies, expectation values, particle distributions, and save wave functions |
| Transition dipole codes | `RG_0S-1P`, `RG_1P-2D`, `RG_1P-2P` | Calculation of off-diagonal matrix elements of the transition dipole moment operator in the length and velocity gauges. |
| Fine structure coupling codes | `RG_0S-2D`, `RG_0S-2P`, `RG_1P-1P`, `RG_2D-2D`, `RG_2P-2D`, `RG_2P-2P` | Calculations of off-diagonal matrix elements of the spin–orbit interaction and non-contact spin–spin interaction |

Each code is a self-contained project. Its directory contains a `Makefile` and a subdirectory `src` with the actual source. The source structure and the structure of the makefiles for the cades within each group are very similar. Each directory with a code also contains an subdirectory `.vscode` with standard JSON configuration files for Microsoft Visual Studio Code (VS Code). These are user-independent and provide the capability to build and debug each code in its directory using VS Code.

Most code directories also contain a subdirectory called `sample_input`. These contain simple examples of input files as well as instructions (in files `README.md`) on how to execute those sample calculations and any additional files or scripts that may be relevant. For the energy and transition dipole codes each sample case directory holds its own `README.md`; for the fine structure coupling codes the instructions are instead placed inside the case's `initial_state_*`, `final_state_*`, and `matelem` subdirectories.

## Root directory files

There are several files located in the project's root directory `ecgpack/`. Their description is provided below:

| File | Description |
| --- | --- |
| `README.md` | This repository manual file. |
| `AGENTS.md` | Configuration rules for AI coding agents and autonomous workflows (e.g. OpenAI Codex). |
| `AUTHORS.md` | List of ECGPACK contributors. |
| `CITATION.cff` | Machine-readable citation metadata in the Citation File Format. GitHub uses it to generate the "Cite this repository" link. |
| `CLAUDE.md` | Configuration and system instructions for Anthropic's Claude Code CLI. |
| `LICENSE.md` | The license under which ECGPACK is distributed (BSD 3-Clause). |
| `THIRD-PARTY-NOTICES.md` | Copyright and license notices for the third-party source code bundled with the codes (the netlib reference BLAS and LAPACK, the PORT optimization routines, and the SLATEC machine constants). |
| `build.bash` | A Bash-script for batch compilation of multiple codes corresponding to a different number of particles, toolchains, configurations, precision, etc. It is convenient for building a large number of different binaries that are later used in production calculations. The generated binaries are automatically moved to directory `/ecgpack/bin`. For more information run this script in a terminal without arguments or read its header. |
| `.code-workspace` | A JSON configuration file for Microsoft Visual Studio Code (VS Code) that contains information used to group separate code project directories into a single, unified workspace that can be opened in VS Code. |

## Compilation and execution

Each code is compiled by going to its directory and running `make` with the appropriate arguments (compiler/toolchain, configuration, working precision, and linear algebra library). More conveniently, one can use the `build.bash` script in the root directory to batch-compile many code variants in one step. The codes can be built with double (fp64), extended (fp80), or quadruple (fp128) precision, and the energy codes can be linked against either the bundled netlib reference BLAS/LAPACK or an optimized library (MKL, OpenBLAS, AOCL, etc.). It is important to note that number of particles is compiled in rather than supplied at runtime. Each binary is an MPI program that is launched in the usual way with `mpirun -np <NPROCS> <BINARYFILE>` from the work directory containing the required input file(s). For full details on the build arguments, precision and performance trade-offs, linear algebra options, the `build.bash` script, binary naming, and execution, see [Compilation and execution](doc/compilation_and_execution.md) in the documentation folder.

## Input file format

The energy codes read and write a single input/output file named `inout.txt` located in the work directory where the code is executed. This file has a certain format and consists of four sections: a header (defining the quantum system and solver parameters), a command list (the sequence of actions to perform), a history (energies obtained at each basis size), and the basis functions themselves. For a full description on the input file format see [Input file format](doc/input_file_format.md) in the documentation folder.

## Use of Microsoft Visual Studio Code

For editing, compiling, debugging, and browsing the code locally, we recommend using Microsoft Visual Studio Code. This ECGPACK repository includes pre-configured JSON settings for the editor. If you are new to VS Code, a quick-start guide is available in the documentation folder: [Use of Microsoft Visual Studio Code](doc/use_of_visual_studio_code.md).

## License [![License: BSD 3-Clause](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](LICENSE.md)

ECGPACK is distributed under the BSD 3-Clause License. The full text is in [`LICENSE.md`](LICENSE.md).

Some amount of source code written by others is bundled with the codes and remains under the terms of its original authors: the netlib reference BLAS and LAPACK (`src/BLAS.f` and `src/LAPACK.f`), the PORT unconstrained minimization routines of David M. Gay (`src/dmng.f`), and the SLATEC machine constants `I1MACH` and `D1MACH` (`src/X1MACH.f90`). The corresponding copyright and license notices are collected in [`THIRD-PARTY-NOTICES.md`](THIRD-PARTY-NOTICES.md).
