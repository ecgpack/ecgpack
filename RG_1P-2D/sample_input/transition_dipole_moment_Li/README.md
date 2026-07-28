# $P^o \rightarrow D^e$ Transition Dipole Moment Calculation (Lithium atom)

## Description

This calculation evaluates the **electric-dipole transition matrix element** between the initial $2~^2P^o$ state and the final $3~^2D^e$ state of the lithium atom using the ECG transition dipole moment code.

The `RG_1P-2D` code calculates only the **transition matrix elements** (the length-form coordinate matrix element and the velocity-form momentum matrix element). 
It does **not** calculate the line strength or oscillator strength directly.

The line strength and oscillator strength are calculated afterward in a separate post-processing step, for example in an Excel worksheet, using the calculated transition matrix elements, transition energy, and the required angular-momentum factors. 
For the definitions of the reduced transition matrix elements, line strengths, and length- and velocity-form oscillator strengths, see:

S. Nasiri, L. Adamowicz, and S. Bubin, *Physical Review A* **112**, 062809 (2025), https://doi.org/10.1103/qrtf-56np.

Only the wave-function files of the initial ($P^o$) and final ($D^e$) states are required when the transition code is executed.
No input file such as `inout.txt` is required at that stage.

The initial and final wave functions must correspond to the **same isotope** and must use the same particle masses and particle ordering.
Otherwise, the transition dipole moment calculation will fail or the resulting matrix element will not be physically meaningful.


## 1. Required files

Two independently optimized ECG wave functions are required:

| State | Role in the transition |
|---|---|
| $P^o$ state | Initial state |
| $D^e$ state | Final state |

The wave functions must be calculated for the same isotope and with consistent physical parameters.


## 2. Save each wave function

In the example directories in the `RG_1P` and `RG_2D` codes, there is a `store_wavefunction` folder. 
The input files in those directories can be used to compute the wave function.
Run the $P^o$-state and $D^e$-state calculations separately and save both wave functions.


## 3. Rename the wave-function files

The `RG_1P-2D` transition code expects the following file names:

| State | Required file name |
|---|---|
| Initial $P^o$ state | `wf_state1.txt` |
| Final $D^e$ state | `wf_state2.txt` |

Rename or copy the saved wave functions accordingly.

Copy the generated wave functions from the initial- and final-state folders into the `transition_dipole_moment_Li` working directory and rename them accordingly.


## 4. Prepare the working directory

The `transition_dipole_moment_Li` folder is the working directory for the transition calculation.
Place the following files in this folder:

```text
RG_1P-2D
wf_state1.txt
wf_state2.txt
```

An `inout.txt` file is not required in this directory because the transition code reads the two wave-function files directly.

After the calculation, the result is presented in the terminal, or can be written to a file (see Section 5).


## 5. Run the transition calculation

Run the `RG_1P-2D` executable from the `transition_dipole_moment_Li` working directory.

For example:

```bash
./RG_1P-2D > out.txt
```

The program reads:

```text
wf_state1.txt
wf_state2.txt
```

and calculates the transition matrix elements for the $P^o\rightarrow D^e$ transition.


## 6. Use the calculated transition matrix elements

The results obtained from `RG_1P-2D` are the transition matrix elements associated with the selected components of the electric-dipole operator: the $z$ component in the length gauge and the $P_z$ component in the velocity gauge.

The orbital reduced transition matrix element is obtained from the calculated component using the appropriate angular-momentum factor.
The line strength is then calculated from the squared magnitude of the reduced matrix element.
Finally, the oscillator strength is calculated using the line strength, transition energy, and initial-state statistical weight.

For example, these calculations may be carried out in an Excel worksheet using:

- the transition matrix elements obtained from the `RG_1P-2D` code;
- the transition energy;
- the statistical weight of the initial state;
- the required $3j$-symbol factor;
- the required $6j$-symbol factor; and
- the length- or velocity-gauge oscillator-strength expression.


## Calculation workflow

```text

Save the initial P-state wave function
             +
Save the final D-state wave function
             ↓
Copy and rename the initial P-state wave function as wf_state1.txt
             ↓
Copy and rename the final D-state wave function as wf_state2.txt
             ↓
Run RG_1P-2D in the transition_dipole_moment_Li folder
             ↓
Generate the resultant output file and obtain the transition matrix elements
             ↓
Calculate the line strength and oscillator strength separately, for example in Excel
```


## Final results

Atomic units are used throughout.

For this sample calculation, the transition energy is

$$
\Delta E =  0.074\,525\,543.
$$

The calculated length-form coordinate matrix element is

$$
\left\langle P^o \middle| z \middle| D^e \right\rangle
=2.219\,088\,687.
$$

The calculated velocity-form momentum matrix element is purely imaginary and is given by

$$
\left\langle P^o \middle| P_z \middle| D^e \right\rangle
=0.169\,427\,277\;i.
$$

The corresponding line strengths and oscillator strengths are evaluated separately during the post-processing step, for example using an Excel worksheet.

The length-form line strength is

$$
S^{\mathrm{L}} = 73.865\,318\,98,
$$

and the velocity-form line strength is

$$
S^{\mathrm{V}} = 0.430\,584\,032.
$$

The corresponding oscillator strengths are

$$
f^{\mathrm{L}} = 0.611\,650\,338,
$$

and

$$
f^{\mathrm{V}} = 0.641\,963\,386.
$$

For comparison, a high-accuracy nonrelativistic reference value for the infinite-nuclear-mass lithium atom, obtained using Hylleraas-type wave functions, is:

$$
f^{\mathrm{L}}_{\mathrm{ref}} = 0.638\,568\,129.
$$


D. Sun, L. Wang, Z.-C. Yan, *Atomic Data and Nuclear Data Tables* **169**, 101793 (2026), https://doi.org/10.1016/j.adt.2026.101793.


The disagreement between the present length- and velocity-form oscillator strengths, as well as their deviations from the reference value, is primarily attributable to incomplete basis-set convergence because only 100 ECG basis functions are used for each state in this sample calculation.