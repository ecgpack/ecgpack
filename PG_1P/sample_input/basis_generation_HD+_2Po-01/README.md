# Description

Sample file `inout.txt` in this directory provides input for generating a very small ECG basis for the lowest $L = 1$ rovibrational state $(v = 0)$ of the $\mathrm{HD}^+$ ion (deuteron, proton and electron, all treated non-adiabatically; the state is the $^2P^o$ counterpart of the $^2S^e$ ground state). The basis functions are $z_i\, r_1^{m}\, \exp[-\mathbf{r}'(L_k L_k' \otimes I_3)\mathbf{r}]$ with the angular prefactor on either the internuclear vector ($z_1$) or the electron ($z_2$); the index $i$ and the even power $m$ of every new function are selected by an exhaustive scan. It starts from 0 basis functions and stops when the basis reaches 10 basis functions. When the basis size equals 5 and 10 functions, the program will save the input in separate files  
`inout_HD+_2Po-01-00005.txt`  
`inout_HD+_2Po-01-00010.txt`  
The calculation with this input file should take about 2 seconds on a single CPU core when double precision is used (gfortran compiler / Intel Core i9-7900X class CPU).
At the end of the calculation the energy should be about -0.596 (a test run with `ECG_RND_SEED=12345` gave -0.5917901 with 5 functions and -0.5963517 with 10 functions; the converged value for this state is -0.597698128 hartree, J. Phys. B 48, 195101 (2015)). With such a small basis every selected function carries the $z_1$ prefactor; the $z_2$ functions enter at the $10^{-8}$ hartree level in large bases.

Note that the end value of the energy may differ slightly from one execution to another because each execution involves stochastic selection of the basis functions. Setting the environment variable `ECG_RND_SEED` to an integer makes the selection reproducible for the same number of MPI processes.

The warning WC0002 about `CURRENT_ENERGY` (1.0E30) and the message `*** WARNING: the previous step left a non-physical energy ***` printed before the first step are expected: 1.0E30 is the placeholder value of a calculation that starts from an empty basis.
