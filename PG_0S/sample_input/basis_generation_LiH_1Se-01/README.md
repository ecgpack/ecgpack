# Description

Sample file `inout.txt` in this directory provides input for generating a very small ECG basis for the ground state of the LiH molecule (lowest $^1 S^e$ state; lithium, proton and four electrons, all treated non-adiabatically). The basis functions have the form $\phi_k=r_1^{2 m_k} \exp [ \mathbf{r}' (A_k \otimes \mathbf{I}) \mathbf{r}]$, where $r_1$ is the internuclear distance. It starts from 0 basis functions and stops when the basis reaches 10 basis functions. When the basis size equals 5 and 10 functions, the program will save the input in separate files  
`inout_LiH_1Se-01-00005.txt`  
`inout_LiH_1Se-01-00010.txt`  
The calculation with this input file should take about 20 seconds on a single CPU core when double precision is used (gfortran compiler / Intel Core i9-7900X).
At the end of the calculation the energy should be about -8.00 (five test runs gave -7.9320 to -7.9456 with 5 functions and -7.9911 to -8.0021 with 10 functions). Such a small basis reproduces only the first decimal figure of the exact non-Born-Oppenheimer energy of LiH, -8.06643920.

Note that the end value of the energy may differ slightly from one execution to another because each execution involves stochastic selection of the basis functions. Setting the environment variable `ECG_RND_SEED` to an integer makes the selection reproducible.

The warning WC0002 about `CURRENT_ENERGY` (1.0E30) is expected: it is the placeholder value of a calculation that starts from an empty basis.
