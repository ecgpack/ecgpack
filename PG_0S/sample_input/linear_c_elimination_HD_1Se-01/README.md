# Removing functions with small linear coefficients

## Description

The sample `inout.txt` file in this directory defines a basis of 10 ECG functions for the ground state of the HD molecule (lowest $^1S^e$ state; deuteron, proton and two electrons, all treated non-adiabatically).
Running it removes every basis function whose linear coefficient is smaller in magnitude than the threshold given in the input, and writes the reduced basis, in the format of `inout.txt`, to a separate file whose name is specified in the input (`basis_elim_lcfn.txt`). `inout.txt` itself is not changed.

With the threshold 0.05, functions 1, 5, and 10 are removed, so the basis shrinks from 10 to 7 functions and the energy rises from -1.1464025232 to -1.1199623315.

On a single CPU core, the calculation takes about 1 second in double precision (gfortran compiler, Intel Core i9-7900X).

## The `ELIM_LCFN` instruction

To reduce the basis, the following instruction is added to `inout.txt`:

```text
ELIM_LCFN G 10 0.05 basis_elim_lcfn.txt
```

The relevant arguments are:

| Argument | Description |
|---|---|
| `ELIM_LCFN` | Keyword for removing functions with small linear coefficients |
| `G` | Eigensolver type (`G` or `Q`; with `I` the step is skipped) |
| `10` | Number of ECG basis functions (must equal the current basis size) |
| `0.05` | Threshold for the magnitude of the linear coefficients of the normalized functions |
| `basis_elim_lcfn.txt` | Name of the file for the reduced basis |

The program stops after writing the file, so this must be the last instruction in `inout.txt`. The stop is made through `MPI_Abort`, so Open MPI prints an `MPI_ABORT was invoked ... with errorcode 1` notice at the end; this is expected.
