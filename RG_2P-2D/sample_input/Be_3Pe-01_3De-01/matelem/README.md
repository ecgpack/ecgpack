# Description

## Input files

The sample files `wf_state0.txt` and `wf_state1.txt` contain the wave functions of the lowest 
beryllium $^3P^{e}$ and $^3D^{e}$ states, respectively.

These files are extracted from the general-format input files using the `RG_2P` and `RG_2D`  codes.

The program requires an input file named `inp.txt` containing a single positional argument:

**Transition type**

* `3P_1D` - offdiagonal matrix element


$$
\langle {}^{3}P_{J=2} | \mathcal{H} | {}^{1}D_{J=2} \rangle.
$$


* `3P_3D` - offdiagonal matrix element

$$
\langle {}^{3}P_{J=1} | \mathcal{H} | {}^{3}D_{J=1} \rangle.
$$

* `1P_3D` - offdiagonal matrix element

$$
\langle {}^{1}P_{J=1} | \mathcal{H} | {}^{3}D_{J=1} \rangle.
$$

* `4P_2D` - offdiagonal matrix element

$$
\langle {}^{4}P_{J=5/2} | \mathcal{H} | {}^{2}D_{J=5/2} \rangle.
$$

Here, $\mathcal{H}$ denotes the spin-dependent part of the Breit–Pauli Hamiltonian. 

### Example

For the given sample files, the following line in `inp.txt` requests the calculation of the spin-dependent contribution to

$$
\langle {}^{3}P_{J=1} | \mathcal{H} | {}^{3}D_{J=1} \rangle.
$$

```text
3P_3D
```

## Execution

Upon execution, the program reads `wf_state0.txt` and `wf_state1.txt` and evaluates the requested matrix elements of the Breit–Pauli Hamiltonian between the corresponding states.

For the supplied sample wave functions, the calculation completes within a few seconds on a single CPU core.

## Output files

The results are written to file `expvals.txt`. An auxiliary file, `spinData.txt`, contains the spin wave functions of the initial and final states.

## Additional notes

1. For calculations of

$$
\langle {}^{3}P_{J=2} | \mathcal{H} | {}^{1}D_{J=2} \rangle,
$$

the input files must contain:

* `wf_state0.txt` — the wave function $|{}^{3}P\rangle$;
* `wf_state1.txt` — the wave function $|{}^{1}D\rangle$.


2. For calculations of 

$$
\langle {}^{3}P_{J=1} | \mathcal{H} | {}^{3}D_{J=1} \rangle.
$$

the input files must contain:

* `wf_state0.txt` — the wave function $|{}^{3}P\rangle$;
* `wf_state1.txt` — the wave function $|{}^{3}D\rangle$.

3. For calculations of 

$$
\langle {}^{1}P_{J=1} | \mathcal{H} | {}^{3}D_{J=1} \rangle.
$$

the input files must contain:

* `wf_state0.txt` — the wave function $|{}^{1}P\rangle$;
* `wf_state1.txt` — the wave function $|{}^{3}D\rangle$.


4. For calculations of 

$$
\langle {}^{4}P_{J=5/2} | \mathcal{H} | {}^{2}D_{J=5/2} \rangle.
$$

the input files must contain:

* `wf_state0.txt` — the wave function $|{}^{4}P\rangle$;
* `wf_state1.txt` — the wave function $|{}^{2}D\rangle$.


The order of these files is important and should not be changed.
