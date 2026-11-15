---
name: run-sample-case
description: Build an ECGPACK energy code if needed and run one of its sample_input cases in a scratch copy under jobs/, then compare the result with the reference values in the case README. Use when asked to validate a code change, reproduce a sample calculation, or check that a build works.
---

# Run a sample case

ECGPACK has no automated test suite. The standard check after a change that
could affect numerical results is to run a shipped example and compare it with
the reference values in that example's `README.md`.

## Steps

1. Choose a case: `ls <CODE>/sample_input/`. Read its `README.md` for the
   expected runtime and reference values. The `basis_generation_*` cases take
   a few seconds.

2. From the repository root, run the helper script of this skill:

   ```bash
   bash .agents/skills/run-sample-case/scripts/run_sample_case.sh <CODE> <case> [nprocs]
   ```

   For example:

   ```bash
   bash .agents/skills/run-sample-case/scripts/run_sample_case.sh RG_0S basis_generation_He_1Se-01
   ```

   The script reads `PARTICLES` from the case's `inout.txt` (the number of
   particles is compiled into each binary), builds the binary with
   `build.bash` if it is missing, copies the case into a new directory under
   `jobs/`, runs it there with `mpirun`, and prints the resulting header
   values and the job directory. It never runs inside `sample_input/`, because
   the codes rewrite `inout.txt` in place.

3. Compare the printed `CURRENT_ENERGY`, or the expectation values in
   `output.log`, with the case `README.md`. State which numbers you compared
   and whether they agree.

4. Delete the job directory when it is no longer needed (`jobs/` is
   gitignored).

## Options

Set these environment variables to change the defaults:

| Variable | Default | Meaning |
| --- | --- | --- |
| `REBUILD` | `0` | `1` forces a rebuild. **Use it after editing source files**, otherwise an existing, stale binary is reused. |
| `MACHINE` | `linux-generic` | `machine` argument of `build.bash` |
| `TOOLCHAIN` | `systemdefault` | `toolchain` argument of `build.bash` |
| `CONFIG` | `release` | `release` or `debug` |
| `PRECISION` | `8` | `8`, `10`, or `16` |
| `LINALG` | `netlib` | linear algebra library |

## Notes

- Cases whose `inout.txt` has a `RANDOM_SEED` line give the same result in
  repeated runs with the same number of MPI processes, apart from rounding
  differences. To compare results before and after a code change, use the
  same `nprocs` in both runs.
- The script handles the energy codes that `build.bash` can build: `RG_0S`,
  `RG_1P`, `RG_2D`, `RG_2P`, and `PG_0S`. `build.bash` does not know `PG_1P`
  yet: build that binary through `PG_1P/Makefile` first, and the script will
  then use it.
- The off-diagonal codes (`RG_0S-1P` and the like) read wave-function files
  instead of `inout.txt`. Follow the `README.md` files in their sample
  directories.
