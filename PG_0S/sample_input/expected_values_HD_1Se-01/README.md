# Description

Sample file `inout.txt` in this directory contains a basis of 10 ECG functions for the lowest $^1 S^e$ state of the HD molecule (deuteron, proton and two electrons, all treated non-adiabatically) and sets up calculations of expectation values for various quantities.
The calculations should take about 2 second on a single CPU core when double precision is used (gfortran compiler / Intel Core i9-7900X).
The expectation values are printed and written to `expvals.txt`; the energy should be -1.1464025232 (a test run gave -1.1464025232015291).
