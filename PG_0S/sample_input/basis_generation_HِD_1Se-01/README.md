# Description

Sample file `inout.txt` in this directory provides input for generating a very small ECG basis for the lowest $^1 S^e$ state of the HD molecule (deuteron, proton and two electrons, all treated non-adiabatically). It starts from 0 basis functions and stops when the basis reaches 10 basis functions. When the basis size equals 5 and 10 functions, the program will save the input in separate files  
`inout_HD_1Se-01-00005.txt`  
`inout_HD_1Se-01-00010.txt`  
The calculation with this input file should take about 10 seconds on a single CPU core when double precision is used (gfortran compiler / Intel Core i9-7900X).
At the end of the calculation the energy should be about -1.15 (a test run gave -1.125353855 with 5 functions and -1.150821504 with 10 functions).

Note that the end value of the energy may differ slightly from one execution to another because each execution involves stochastic selection of the basis functions. Setting the environment variable `ECG_RND_SEED` to an integer makes the selection reproducible.

The warning WC0002 about `CURRENT_ENERGY` (1.0E30) is expected: it is the placeholder value of a calculation that starts from an empty basis.
