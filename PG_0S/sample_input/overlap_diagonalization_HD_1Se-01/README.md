# Spectrum of the overlap matrix

## Description

The sample `inout.txt` file in this directory defines a basis of 10 ECG functions for the ground state of the HD molecule (lowest $^1S^e$ state; deuteron, proton and two electrons, all treated non-adiabatically).
Running it diagonalizes the overlap matrix of the normalized basis functions (the Hamiltonian is not involved), prints its smallest and largest eigenvalues and its condition number, and writes the complete spectrum to a separate file whose name is specified in the input (`overlap.txt`). The basis is not changed.

A smallest eigenvalue close to zero means that the basis is nearly linearly dependent; the condition number shows how many decimal digits this costs. For this basis the smallest eigenvalue should be 1.3272328489E-02, the largest 4.8203527431, and the condition number 363.19, so the basis is well conditioned.

On a single CPU core, the calculation takes less than 1 second in double precision (gfortran compiler, Intel Core i9-7900X).

## The `OVERLAP_D` instruction

To obtain the spectrum, the following instruction is added to `inout.txt`:

```text
OVERLAP_D G 10 overlap.txt
```

The relevant arguments are:

| Argument | Description |
|---|---|
| `OVERLAP_D` | Keyword for diagonalizing the overlap matrix |
| `G` | Eigensolver type (only `G` is available; with `I` or `Q` the step is skipped) |
| `10` | Number of ECG basis functions |
| `overlap.txt` | Name of the file for the complete spectrum (optional; `overlap.txt` if omitted) |

The basis-set size given in the `OVERLAP_D` instruction must match the current basis size.
