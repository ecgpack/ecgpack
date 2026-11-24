# Generator and unit test of `qrlinalg.f90`

`src/qrlinalg.f90` is the QR layer of eigensolver method `Q` (see
`../Q_method_design.md`, sections 3 and 4). The file is byte-identical in the
seven energy codes. It holds two state types with the same procedures:

- `qr_real_state` (procedures `real_*`), used by the six real codes; and
- `qr_complex_state` (procedures `complex_*`), the complex Hermitian twin, used
  by `CG_0S`.

**The complex half is generated from the real half.** This directory holds the
generator and the standalone unit test of the module, so that the module can be
changed and checked without any other material.

| File | Purpose |
| --- | --- |
| `gen_complex.py` | The generator. Also recovers the real source from an installed module and checks that a module is up to date. |
| `gen_complex_text.py` | Text that the generator inserts into the declaration part (module-header paragraph, constant, timing model, `ZGEQRF`/`ZUNGQR` interfaces, banner) and the edits of the state type. Executed by `gen_complex.py`; not a program. |
| `test_qrlinalg.f90` | Standalone MPI unit test and timing driver of both state types. |
| `unit_test.sh` | Builds the module and the driver and runs the unit suite at several process counts. |

The scripts need Python 3 (standard library only) and, for the unit test,
gfortran with an MPI compiler wrapper.

## The rule

Edit the real half, regenerate, never edit the complex half by hand. Then copy
the regenerated file to all seven codes.

"The real half" is everything in `qrlinalg.f90` that is not listed under
"What is generated" below: the module header, the declaration part, type
`qr_real_state`, the type-independent `qr_*` procedures, and the `real_*`
procedures.

## What is generated

A generated module is its real source plus pure insertions. No line of the
real source is changed or removed. The insertions are, in file order:

1. a paragraph of the module header ("Complex Hermitian twin. ..."), directly
   before the `module qrlinalg` statement;
2. the constant `QR_DIAGONAL_TOLERANCE` and the timing model
   `complex_factor_model` with their comments;
3. the statements `public :: qr_complex_state` and
   `public :: complex_factor_model`;
4. the interfaces of `ZGEQRF` and `ZUNGQR`;
5. type `qr_complex_state`, derived from `qr_real_state`;
6. after the last real procedure: a banner, the `complex_*` procedures, each
   derived from its `real_*` counterpart, `complex_dotu` (a second procedure
   derived from `real_dot`), and a few procedures that have no real
   counterpart (`EXTRA_PROCEDURES` in `gen_complex.py`).

A complex procedure is produced from the real one in three steps:

1. procedure-specific substitutions (`EDITS` and `COMMENT_EDITS` in
   `gen_complex.py`: declarations, conjugate transposes, Hermitian fills, MPI
   counts, comments that the generic rules cannot derive). Each substitution
   states how many times it must match. If it matches a different number of
   times the generator stops and names the procedure and the text. A change of
   the real half that touches such a place is therefore noticed and is not
   mistranslated silently;
2. generic renaming (`real_*` to `complex_*`, `DGEQRF` to `ZGEQRF`, and so on);
3. generic rewording of comments ("symmetric" to "Hermitian", `^T` to `^H`, and
   so on), followed by rewrapping of the lines that became too long.

## Changing the module

Run the commands from the repository root.

```bash
G=doc/devnotes/qrlinalg_generator

# 1. Before editing: the installed file must be what the generator produces.
python3 $G/gen_complex.py --check RG_0S/src/qrlinalg.f90

# 2. Edit the real half of RG_0S/src/qrlinalg.f90.

# 3. Regenerate the complex half in place and review the change.
python3 $G/gen_complex.py --regenerate RG_0S/src/qrlinalg.f90 RG_0S/src/qrlinalg.f90
git diff RG_0S/src/qrlinalg.f90

# 4. Run the unit test (next section).

# 5. Install the file in the other six codes and confirm that all are identical.
for c in RG_1P RG_2D RG_2P PG_0S PG_1P CG_0S; do cp RG_0S/src/qrlinalg.f90 $c/src/; done
md5sum */src/qrlinalg.f90
```

`--regenerate MODULE OUTPUT` extracts the real source contained in `MODULE`,
generates the full module from it, and writes it to `OUTPUT`. The complex half
of `MODULE` is ignored, so it may be out of date. `--check MODULE` does the
same in memory and compares the result with `MODULE`: exit status 0 when the
file is byte for byte what the generator produces from its own real half, and
1 otherwise.

When the generator stops with `expected N occurrence(s) of ..., found M`, the
edit of the real half changed a line that a procedure-specific substitution
refers to. Adapt that substitution in `gen_complex.py` (or add one for a new
procedure whose complex version needs more than renaming) and regenerate. A new
real procedure that needs no substitution is translated by the generic steps
alone. A procedure whose name starts with `qr_` is type independent and is not
duplicated.

The generator prints the numbers of lines that are longer than 96 characters,
if any. Shorten such a line in the real half or in the substitution that
produced it.

### Working with a separate real source

The module was developed with the real half in a file of its own. That form is
still available:

```bash
python3 $G/gen_complex.py --extract RG_0S/src/qrlinalg.f90 /tmp/qrlinalg_real.f90
python3 $G/gen_complex.py /tmp/qrlinalg_real.f90 /tmp/qrlinalg.f90
diff /tmp/qrlinalg_real.f90 /tmp/qrlinalg.f90 | grep -c '^<'    # prints 0: pure insertions
cmp /tmp/qrlinalg.f90 RG_0S/src/qrlinalg.f90                    # no output: identical
```

This is also the round-trip check of the generator itself: extraction followed
by generation must reproduce the installed file. The real source is not kept
in the repository, because it is contained in the installed file.

### Changing the inserted text

The extraction finds insertions 1, 2, and 4 by their exact text in
`gen_complex_text.py`, insertion 3 by its statements, and insertions 5 and 6 by
their first and last lines. To change the text of insertion 1, 2, or 4, or the
banner of insertion 6, first extract the real source with the unchanged
scripts, then edit `gen_complex_text.py`, then generate from the extracted
file (the two commands above) and copy the result over
`RG_0S/src/qrlinalg.f90`.

## Unit test

`test_qrlinalg.f90` builds replicated matrices H and S, drives each state type
through random sequences of `begin_active`, `replace_symmetric`,
`append_symmetric`, `truncate`, and flushes, and compares every bordered solve
with a dense reference and the eigenpairs with a dense inverse iteration. It
also tests both fresh-factorization algorithms directly, the shifted quotient
against an evaluation in a wider real kind, the error paths, and the
conventions of the complex state. Its header comment lists the command-line
arguments, including the `timing` and `factortime` modes.

```bash
doc/devnotes/qrlinalg_generator/unit_test.sh 8 10                     # release, 1, 2, 3 processes
NPROCS="1 2 3 5 14" doc/devnotes/qrlinalg_generator/unit_test.sh 8 10
NPROCS=1 doc/devnotes/qrlinalg_generator/unit_test.sh 16              # quadruple precision: one process
MODE=debug NPROCS="1 3" doc/devnotes/qrlinalg_generator/unit_test.sh 8 10
TEST_ARGS=factor=serial      doc/devnotes/qrlinalg_generator/unit_test.sh 8
TEST_ARGS=factor=distributed doc/devnotes/qrlinalg_generator/unit_test.sh 8
```

The script compiles `qrlinalg.f90`, `wp_def_<PREC>.f90`, `X1MACH.f90`,
`BLAS.f`, and `LAPACK.f` from `RG_0S/src` (another code with `CODE_SRC=...`)
with `-DQRLINALG_TESTING`, which makes the components of the state types
public for the white-box tests. That symbol must stay undefined in an ECGPACK
build. Nothing is written into the source tree: objects, the executable, logs,
and traces go to `jobs/qrlinalg_unit/p<PREC>[_dbg]/` (another place with
`BUILD_ROOT=...`). The header of the script lists all settings.
`TEST_ARGS=quick` shortens the random sequences and omits the largest cases.

Every run must end with `RESULT: ALL TESTS PASSED`; the script prints one line
per run and returns a nonzero status otherwise. The first process count of
`NPROCS` writes a trace of the checked values, and every other count is
compared with it within the tolerances stored in the traces.

A change of the module is complete when the unit test passes

- at 1, 2, 3, 5, and 14 processes in precisions 8 and 10 (release);
- at one process in precision 16;
- in strict debug builds (`MODE=debug`) of precisions 8 and 10 at 1 and 3
  processes; and
- each of these with the automatic choice of the factorization algorithm, with
  `factor=serial`, and with `factor=distributed`.

Then validate the energy codes as section 15 of `../Q_method_design.md`
requires.
