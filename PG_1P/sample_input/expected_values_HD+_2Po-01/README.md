# Description

Sample file `inout.txt` in this directory contains a basis of 10 ECG functions for the lowest $L = 1$ rovibrational state of the $\mathrm{HD}^+$ ion (deuteron, proton and electron, all treated non-adiabatically) and sets up the calculation of expectation values.
The calculation should take about 1 second on a single CPU core when double precision is used (gfortran compiler / Intel Core i9-7900X class CPU).
The energy, the norm of the eigenvector and the interparticle distances $\langle r_{ij}\rangle$ and $\langle r_{ij}^2\rangle$ are printed and written to `expvals.txt`; the energy should be -0.5955901654 (a test run gave -0.59559016544462284, $\langle r_1\rangle$ = 2.0429556229 for the deuteron-proton distance, $\langle r_2\rangle$ = 1.6794627468 for the deuteron-electron distance and $\langle r_{12}\rangle$ = 1.6817790047 for the proton-electron distance).

The basis contains only $z_1$ functions, for which the distance formulas are those of the $z_1 r_1^m$ basis of the original code; functions with a $z_i$ prefactor, $i > 1$, are handled with the general formulas of `ExpValuesXi` (module matelem).
