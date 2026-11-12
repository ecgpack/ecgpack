# Description

Sample file `inout.txt` in this directory provides input for generating an ECG basis for the lowest $^2 D^e$ state of Li atom (Li-7 isotope). It starts from 0 basis functions and stops when the basis reaches 100 basis functions. When the basis size equals 5, 10, 20, 50, and 100 functions, the program will save the input in separate files named  
`inout_Li_2De-01-00005.txt`  
`inout_Li_2De-01-00010.txt`  
`inout_Li_2De-01-00020.txt`  
`inout_Li_2De-01-00050.txt`  
`inout_Li_2De-01-00100.txt`  
The calculation with this input file should take roughly 2-3 minutes on a single CPU core when double precision is used (gfortran compiler / AMD Ryzen 7 7800X3D).  
At the end of the calculation the energy should reproduce four decimal figures in the exact value of -7.334927306.

The basis functions are selected stochastically. The `RANDOM_SEED` line in `inout.txt` fixes the seed of the random number generators, and the codes use their own generators, whose sequences do not depend on the compiler. Repeated runs with the same number of MPI processes therefore reproduce the same energies on any machine and with any compiler, apart from rounding differences, which the stochastic selection and the optimization can amplify (especially in longer runs). To get a different random selection in every run, remove the `RANDOM_SEED` line or replace the seed with `SYSTEM_CLOCK` (see `doc/input_file_format.md`).
