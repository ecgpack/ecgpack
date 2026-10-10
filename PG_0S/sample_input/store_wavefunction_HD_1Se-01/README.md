# Writing the wave function

## Description

The sample `inout.txt` file in this directory defines a basis of 10 ECG functions for the ground state of the HD molecule (lowest $^1S^e$ state; deuteron, proton and two electrons, all treated non-adiabatically).
Running it writes the Hamiltonian and overlap matrices, the eigenvector, and the wave function to separate files whose names are specified in the input. The energy should be -1.1464025232.

On a single CPU core, the calculation takes less than 1 second in double precision (gfortran compiler, Intel Core i9-7900X).

## The `SAVE_HSWF` instruction

To produce the files, the following instruction is added to `inout.txt`:

```text
SAVE_HSWF I 10 H_full.txt S_full.txt eigvec.txt wavefunction.txt
```

The relevant arguments are:

| Argument | Description |
|---|---|
| `SAVE_HSWF` | Keyword for writing the matrices and the wave function |
| `I` | Eigensolver type (`G`, `I`, or `Q`) |
| `10` | Number of ECG basis functions |
| `H_full.txt` | Hamiltonian matrix of the normalized basis functions (one element per line, preceded by its two indices) |
| `S_full.txt` | Overlap matrix of the normalized basis functions (same layout) |
| `eigvec.txt` | Linear coefficients of the normalized basis functions |
| `wavefunction.txt` | Wave function: the system header, then one line per basis function with its linear coefficient, the power $2m_k$ of $r_1$, and its nonlinear parameters |

Any of the four file names may be set to `none` to skip that file.
The basis-set size given in the `SAVE_HSWF` instruction must match the size of the wave function being saved.
