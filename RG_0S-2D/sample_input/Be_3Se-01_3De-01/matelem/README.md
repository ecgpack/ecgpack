# Description

## Input files

The sample files `wf_state0.txt` and `wf_state1.txt` contain the wave functions of the lowest beryllium $^{3}S^{e}$ and $^3D^{e}$ states, respectively.

These files are extracted from the general-format input files using the `RG_0S` and `RG_2D` codes.

The program calculates the spin-dependent part of the Breit-Pauli Hamiltonian which for the given transition includes only non-contact spin-spin interaction:
$$
\langle {}^{3}D_{J=1} | \mathcal{H} | {}^{3}S_{J=1} \rangle.
$$

## Execution

Upon execution, the program reads `wf_state0.txt` and `wf_state1.txt` and evaluates the matrix elements of the Breit–Pauli Hamiltonian between the corresponding states.

For the supplied sample wave functions, the calculation completes within a few seconds on a single CPU core.

## Output files

The results are written to file `expvals.txt`.

An auxiliary file, `spinData.txt`, contains the spin wave functions of the initial and final states used in the calculation.

## Additional notes

The input files must contain:

* `wf_state0.txt` — the wave function $|{}^{3}S\rangle$;
* `wf_state1.txt` — the wave function $|{}^{3}D\rangle$.

The order of these files is important and should not be changed.
