# Description

Sample file `inout.txt` in this directory provides input for generating a very small ECG basis for the ground state of He atom (lowest $^1 S^e$ state of He-4 isotope). It starts from 0 basis functions and stops when the basis reaches 10 basis functions. When the basis size equals 5 and 10 functions, the program will save the input in separate files  
`inout_He_1Se-01-00005.txt`  
`inout_He_1Se-01-00010.txt`  
The calculation with this input file should take about 1 second on a single CPU core when double precision is used (gfortran compiler / AMD Ryzen AI 9 HX 370).
At the end of the calculation the energy should reproduce three decimal figures in the exact value of -2.9033045577.

The basis functions are selected stochastically. The `RANDOM_SEED` line in `inout.txt` fixes the seed of the random number generators, and the codes use their own generators, whose sequences do not depend on the compiler. Repeated runs with the same number of MPI processes therefore reproduce the same energies on any machine and with any compiler, apart from rounding differences, which the stochastic selection and the optimization can amplify (especially in longer runs). To get a different random selection in every run, remove the `RANDOM_SEED` line or replace the seed with `SYSTEM_CLOCK` (see `doc/input_file_format.md`).  
