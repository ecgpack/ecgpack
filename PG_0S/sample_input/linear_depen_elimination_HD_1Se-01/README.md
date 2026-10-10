# Removing linearly dependent functions

## Description

The sample `inout.txt` file in this directory defines a basis of 10 ECG functions for the ground state of the HD molecule (lowest $^1S^e$ state; deuteron, proton and two electrons, all treated non-adiabatically).
Running it looks for pairs of normalized basis functions whose overlap exceeds the threshold given in the input, removes the higher-numbered function of each such pair, and writes the reduced basis, in the format of `inout.txt`, to a separate file whose name is specified in the input (`basis_elim_lnd1.txt`). `inout.txt` itself is not changed.

With the threshold 0.9, functions 2 and 7 (overlap 0.970) form the only such pair, so function 7 is removed: the basis shrinks from 10 to 9 functions, the largest overlap falls from 0.970 to 0.729, and the energy rises from -1.1464025232 to -1.1320519170.

On a single CPU core, the calculation takes about 1 second in double precision (gfortran compiler, Intel Core i9-7900X).

## The `ELIM_LND1` instruction

To reduce the basis, the following instruction is added to `inout.txt`:

```text
ELIM_LND1 G 10 0.9 basis_elim_lnd1.txt
```

The relevant arguments are:

| Argument | Description |
|---|---|
| `ELIM_LND1` | Keyword for removing linearly dependent functions (pair overlaps) |
| `G` | Eigensolver type (`G` or `Q`; with `I` the step is skipped) |
| `10` | Number of ECG basis functions (must equal the current basis size) |
| `0.9` | Threshold for the magnitude of the overlap of two normalized functions |
| `basis_elim_lnd1.txt` | Name of the file for the reduced basis |

The program stops after writing the file, so this must be the last instruction in `inout.txt`. The stop is made through `MPI_Abort`, so Open MPI prints an `MPI_ABORT was invoked ... with errorcode 1` notice at the end; this is expected.
