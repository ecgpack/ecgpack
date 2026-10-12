# Checking the basis health

## Description

The sample `inout.txt` file in this directory defines a basis of 10 ECG functions for the ground state of the HD molecule (lowest $^1S^e$ state; deuteron, proton and two electrons, all treated non-adiabatically).
Running it writes the Hamiltonian and overlap matrices of the unnormalized basis functions to separate files whose names are specified in the input (`H_raw.txt` and `S_raw.txt`), and prints a basis health table, which is also written to `basis_health.txt`. No eigenvalue problem is solved.

The table shows, for every basis function, how much of its self-overlap is lost to cancellation when the symmetry projector (Young operator) acts on it: the cancellation factor $C$, the number of decimal digits lost, and the relative error that remains. In this basis $C=1$ for every function, so no digits are lost. Functions with $C>10^4$ are counted at the end of the table; such functions are rejected during basis generation.

On a single CPU core, the calculation takes less than 1 second in double precision (gfortran compiler, Intel Core i9-7900X).

## The `SAVE_HS_R` instruction

To write the matrices and the table, the following instruction is added to `inout.txt`:

```text
SAVE_HS_R I 10 H_raw.txt S_raw.txt
```

The relevant arguments are:

| Argument | Description |
|---|---|
| `SAVE_HS_R` | Keyword for writing the unnormalized matrices and the basis health table |
| `I` | Eigensolver type (`G`, `I`, or `Q`) |
| `10` | Number of ECG basis functions |
| `H_raw.txt` | Hamiltonian matrix of the unnormalized basis functions (`none` skips it) |
| `S_raw.txt` | Overlap matrix of the unnormalized basis functions (`none` skips it) |

The matrix files have the same layout as those of `SAVE_HSWF` (one element per line, preceded by its two indices). The normalized matrices, which are the ones used in the eigenvalue problem, are written by `SAVE_HSWF` (see `store_wavefunction_HD_1Se-01`).
