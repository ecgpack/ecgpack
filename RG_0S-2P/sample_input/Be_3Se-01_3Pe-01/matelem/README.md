# Description

## Input files

The sample files `wf_state0.txt` and `wf_state1.txt` contain the wave functions of the lowest beryllium $^3S^{e}$ and $^3P^{e}$ states, respectively.

These files are extracted from the general-format input files using the `RG_0S` and `RG_2P` codes.

The program requires an input file named `inp.txt` containing a single positional argument:

**Transition type**

* `3P_1S` - offdiagonal matrix element


$$
\langle {}^{3}P_{J=0} | \mathcal{H} | {}^{1}S_{J=0} \rangle.
$$


* `3P_3S` - offdiagonal matrix element

$$
\langle {}^{3}P_{J=1} | \mathcal{H} | {}^{3}S_{J=1} \rangle.
$$


* `4P_2S` - offdiagonal matrix element

$$
\langle {}^{4}P_{J=1/2} | \mathcal{H} | {}^{2}S_{J=1/2} \rangle.
$$

Here, $\mathcal{H}$ denotes the spin-dependent part of the Breit–Pauli Hamiltonian.

### Example

For the given sample files, the following line in `inp.txt` requests the calculation of the spin-dependent matrix elements:

```text
3P_3S
```

## Execution

Upon execution, the program reads `wf_state0.txt` and `wf_state1.txt` and evaluates the requested matrix elements of the Breit–Pauli Hamiltonian between the corresponding states.

For the supplied sample wave functions, the calculation completes within a few seconds on a single CPU core.

## Output files

The results are written to file `expvals.txt`.
An auxiliary file, `spinData.txt`, contains the spin wave functions of the initial and final states used in the calculation.

## Additional notes

1. For calculations of

$$
\langle {}^{3}P_{J=0} | \mathcal{H} | {}^{1}S_{J=0} \rangle.
$$

the input files must contain:

* `wf_state0.txt` — the wave function $|{}^{1}S\rangle$;
* `wf_state1.txt` — the wave function $|{}^{3}P\rangle$.


2. For calculations of

$$
\langle {}^{3}P_{J=1} | \mathcal{H} | {}^{3}S_{J=1} \rangle.
$$

the input files must contain:

* `wf_state0.txt` — the wave function $|{}^{3}S\rangle$;
* `wf_state1.txt` — the wave function $|{}^{3}P\rangle$.

3. For calculations of

$$
\langle {}^{4}P_{J=1/2} | \mathcal{H} | {}^{2}S_{J=1/2} \rangle.
$$

the input files must contain:

* `wf_state0.txt` — the wave function $|{}^{2}S\rangle$;
* `wf_state1.txt` — the wave function $|{}^{4}P\rangle$.

The order of these files is important and should not be changed.
