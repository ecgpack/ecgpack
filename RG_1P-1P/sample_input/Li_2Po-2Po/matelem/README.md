# Description

## Input files

The sample files `wf_state0.txt` and `wf_state1.txt` contain the wave functions of the lowest and second-lowest lithium $^2P^{o}$ states, respectively.

These files are extracted from the general-format input files using the `RG_1P` code.

The program requires an input file named `inp.txt` containing two positional arguments separated by whitespace:

1. **Hamiltonian component**

* `SPIN` — calculate only the spin-dependent part of the Breit–Pauli Hamiltonian.
* `SCALAR` — calculate only the scalar relativistic part.
* `ALL` — calculate both contributions.

2. **Transition type**

* `XP_XP` — offdiagonal matrix element

$$
\langle {}^{X}P_{J} | \mathcal{H} | {}^{X}P_{J} \rangle.
$$

* `3P_1P` — offdiagonal matrix element

$$
\langle {}^{3}P_{J=1} | \mathcal{H} | {}^{1}P_{J=1} \rangle.
$$

Here, $\mathcal{H}$ denotes either the scalar relativistic or spin-dependent part of the Breit–Pauli Hamiltonian. In the first matrix element, the value of $J$ corresponds to the maximum possible value: $J=S+L$.

### Example

The following `inp.txt` file requests the calculation of the spin-dependent contribution to

$$
\langle {}^{2}P_{J=3/2} | \mathcal{H} | {}^{2}P_{J=3/2} \rangle.
$$

```text
SPIN XP_XP
```

## Execution

Upon execution, the program reads `wf_state0.txt` and `wf_state1.txt` and evaluates the requested matrix elements of the Breit–Pauli Hamiltonian between the corresponding states.

For the supplied sample wave functions, the calculation completes within a few seconds on a single CPU core.

## Output files

The results are written to:

* `expvals_scalar.txt` — scalar relativistic matrix elements;
* `expvals_spin.txt` — spin-dependent matrix elements.

An auxiliary file, `spinData.txt`, contains the spin wave functions of the initial and final states used in the calculation.

## Additional notes

For calculations of

$$
\langle {}^{3}P_{J=1} | \mathcal{H} | {}^{1}P_{J=1} \rangle,
$$

the input files must contain:

* `wf_state0.txt` — the singlet-state wave function $|{}^{1}P\rangle$;
* `wf_state1.txt` — the triplet-state wave function $|{}^{3}P\rangle$.

The order of these files is important and should not be changed.
