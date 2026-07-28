# Using Visual Studio Code with ECGPACK (Linux/WSL)

This is a short guide to editing, building, running, and debugging the ECGPACK codes from
[Visual Studio Code](https://code.visualstudio.com/)
in a Linux or Windows Subsystem for Linux (WSL) environment. It is meant only as a quick start; for the details of the build system and input files see
[compilation_and_execution.md](compilation_and_execution.md) and
[input_file_format.md](input_file_format.md).

Every code directory ships with a `.vscode/` subdirectory containing user-independent `settings.json`, `tasks.json`, and `launch.json` files, so most of the configuration described below works out of the box.

## Prerequisites

You need a working build/run toolchain *before* VS Code can do anything useful. If you do not have a toolchain that includes a Fortran compiler and MPI library you can install it from your Linux distro's official repositories. On a typical Ubuntu/Debian system (native or under WSL) the installation of a minimal set of relevant packages looks as follows:

```bash
sudo apt install build-essential gfortran make gdb openmpi-bin libopenmpi-dev libblas-dev liblapack-dev
pip install fortls     # the Fortran language server (used for code navigation)
```

* **VS Code itself** — install from <https://code.visualstudio.com/>.
* **WSL users** — install VS Code on Windows, then the *WSL* extension (below). Launch VS Code from inside your WSL shell with `code .`; the editor runs on Windows while compilers, `gdb`, and MPI all run inside Linux.
* The build toolchain (`gfortran`, `make`, MPI) and the `gdb` debugger must be the same ones used on the command line — VS Code just drives them.

### Recommended VS Code extensions

Install these from the Extensions view (keyboard shortcut `Ctrl+Shift+X`):

| Extension | ID | Purpose |
| :--- | :--- | :--- |
| [Modern Fortran](https://marketplace.visualstudio.com/items?itemName=fortran-lang.linter-gfortran) | `fortran-lang.linter-gfortran` | Syntax highlighting, linting, formatting, and `fortls` integration |
| [C/C++](https://marketplace.visualstudio.com/items?itemName=ms-vscode.cpptools) | `ms-vscode.cpptools` | Provides the `cppdbg`/GDB debugging back-end used by `launch.json` |
| [WSL](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-wsl) | `ms-vscode-remote.remote-wsl` | Required only when running VS Code on Windows against a WSL distro |

The Modern Fortran extension relies on the `fortls` language server (installed above); the shipped `settings.json` already points `fortran.fortls.path` to it and associates `*.f90` with free-form and `*.f` with fixed-form Fortran.

In addition to the above three required extensions you can install some optional ones that might be useful depending on your exact workflow. Examples include
[Remote - SSH](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-ssh),
[markdownlint](https://marketplace.visualstudio.com/items?itemName=DavidAnson.vscode-markdownlint),
[Makefile Tools](https://marketplace.visualstudio.com/items?itemName=ms-vscode.makefile-tools),
[Claude Code for VS Code](https://marketplace.visualstudio.com/items?itemName=anthropic.claude-code),
[Codex - OpenAI's coding agent](https://marketplace.visualstudio.com/items?itemName=openai.chatgpt).

## Opening the project

You can work with a single code or with the whole collection:

* **One code at a time** — open just that code directory, e.g. `File ▸ Open Folder… ▸ RG_0S/`. This is the simplest setup: its `.vscode/` tasks and launch configurations apply directly.
* **The whole repository** — open the multi-folder workspace file [`.code-workspace`](../.code-workspace) in the repository root (`File ▸ Open Workspace from File…`). It lists every code directory as a separate folder and excludes the `release/`, `debug/`, and `bin/` directories from search and language-server indexing.

Instead of using the `File` menu, you can also open a directory straight from the terminal by typing `code` followed by the directory name. In particular, running

```bash
code .
```

in the ECGPACK root directory opens the whole repository in VS Code, and running it inside a single code directory (e.g. `RG_0S/`) opens just that code. This works both on native Linux and under WSL (with the *WSL* extension installed, as described above). Note that `code .` in the root directory opens the repository as a plain folder; to get the multi-folder layout described above you still need to open the `.code-workspace` file (either through `File ▸ Open Workspace from File…` or by typing `code .code-workspace`).

When opened correctly the Fortran files get syntax highlighting and you can jump to definitions/references across modules (this is `fortls` at work).

## Building

Building is done through VS Code *tasks*, which simply invoke each code's `Makefile` (see [compilation_and_execution.md](compilation_and_execution.md)). Run `Terminal ▸ Run Build Task…` (`Ctrl+Shift+B`) and pick one of the predefined tasks, for example a debug or release build with gfortran. Object and module files are placed in the `debug/` or `release/` subdirectory of the code, and the resulting binary is named after the chosen configuration.

> **Tip.** Always build the **debug** configuration before stepping through code
> in the debugger — it compiles with bounds/uninitialized/FPE checks and full
> debug symbols (`-g`), whereas the release build is optimized (`-O3
> -march=native`) and hard to step through. The debug build is also the one the
> launch configurations expect.

You can edit the commands in `.vscode/tasks.json` if you need a different machine, precision (`PREC`), or linear-algebra back-end (`LINALG`); the arguments map one-to-one onto the `make` arguments documented in [compilation_and_execution.md](compilation_and_execution.md). Cleaning tasks (`make cleandebug` / `cleanrelease` / `cleanest`) are also provided.

## Running and debugging (single MPI process)

Open the **Run and Debug** view (`Ctrl+Shift+D`), choose a configuration from the dropdown at the top (e.g. `run_RG_0S_debug_gfortran_8_netlib`), and press **F5**. Each launch configuration first runs its matching build task and then starts the binary under `gdb`.

A few practical points:

* **Working directory.** The program is launched with its `cwd` set to the code directory, so it reads `inout.txt` (energy codes) or the wave-function files (off-diagonal codes) from there. Put an input file for a small, fast test case in that directory before launching, or edit the `cwd` field in `launch.json` to point at a job directory.
* **Single process only.** The launch configurations start the executable *directly*, not through `mpirun`, so the program runs as a single MPI rank. This is exactly what you want for interactive debugging — VS Code's debugger attaches to one process. (To run in parallel, use `mpirun -np <N> …` from a terminal instead; debugging many ranks from the GUI is not supported here.)
* **Breakpoints.** Click in the gutter to the left of a line number, or press `F9`. Use `F5` (continue), `F10` (step over), `F11` (step into), and `Shift+F11` (step out). The **Debug Console** at the bottom accepts raw `gdb` commands prefixed with `-exec`, e.g. `-exec print Glob_n`.

## Inspecting variables (locals and module-level globals)

While stopped at a breakpoint, the **Variables** panel (top of the Run and Debug view) shows the **Locals** of the current subroutine/function automatically, plus the routine's arguments. Hovering over a variable in the editor shows its value inline.

ECGPACK keeps most of its state in module-level variables (the `Glob_*` names in `globvars.f90`, particle counts and precision in `wp_def_*.f90`, etc.). These are **not** local to the current routine, so they do not appear under *Locals*. To inspect them:

* **Watch panel.** In the **Watch** panel press `+` and type the variable name. Recent `gdb` versions understand Fortran module scoping, so `globvars::Glob_n` (i.e. the format `module::variable`) works. If a bare name such as `Glob_n` is ambiguous or not found, qualify it with the module name.
* **Debug Console.** Equivalently, type in the console:

  ```text
  -exec print globvars::Glob_n
  -exec print Glob_Hkl          # bare name, if unambiguous and in scope
  ```

  If your `gdb` does not accept the `module::variable` form, gfortran exposes the mangled symbol `__<module>_MOD_<variable>` (all lowercase), e.g. `-exec print __globvars_MOD_glob_n`.
* Array slices and individual elements can be printed too, e.g. `-exec print Glob_YMatr(1:3,1:3)`, which is handy for inspecting the Hamiltonian and overlap matrices.

For this to work the variable must actually be live (allocated and in scope at the current execution point) and the code must have been built with the **debug** configuration so that full symbol information is available.
