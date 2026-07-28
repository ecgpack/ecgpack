# Description

Sample file `inout.txt` in this directory provides input for generating a very small ECG basis for the lowest $^1 P^o$ state of He atom (He-4 isotope). It starts from 0 basis functions and stops when the basis reaches 10 basis functions. When the basis size equals 5 and 10 functions, the program will save the input in separate files  
`inout_He_1Po-01-00005.txt`  
`inout_He_1Po-01-00010.txt`  
The calculation with this input file should take about 1 second on a single CPU core when double precision is used (gfortran compiler / AMD Ryzen AI 9 HX 370).
At the end of the calculation the energy should reproduce three decimal figures in the exact value of -2.123545654.

Note that the end value of the energy may differ slightly from one execution to another because each execution involves stochastic selection of the basis functions.  
