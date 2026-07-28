# Compilation and execution

## Compilation

To compile any selected code (e.g. `RG_0S`) using a specific toolchain and options one can go to the code directory (e.g. `RG_0S/`) and run the `make` command with the corresponding arguments. Please look at the header of the `Makefile` in the code directory as it provides a brief explanation of the arguments.

Here are examples of how one can compile ECGPACK codes:

```bash
cd RG_0S
make release COMPILER=gfortran MACHINE=linux-generic PREC=10 LINALG=netlib EXEFILE=mybinaryfile 
```

```bash
cd RG_1P
make debug COMPILER=ifort MACHINE=linux-generic PREC=8 LINALG=mkl EXEFILE=mybinaryfile 
```

```bash
cd RG_2D
make release COMPILER=ifx PREC=16 LINALG=netlib EXEFILE=mybinaryfile 
```

```bash
cd RG_2P
make release COMPILER=nvfortran LINALG=lblas EXEFILE=mybinaryfile 
```

Note that even though the option `COMPILER` here specifies a Fortran compiler name, what is actually called under the hood is the corresponding MPI wrapper (`mpif90`, `mpiifort`, etc).  

To clean the source directory of any object and module files one can run `make clean` or `make cleaner` commands.

### Precision

Each code in the ECGPACK collection can be compiled with a different working precision (real kind) specified by the user (provided there is compiler and hardware architecture support) as an argument in the `make` command:

* `PREC=8` - fp64 (double precision). This is universally supported by all Fortran compilers on pretty much any hardware.
* `PREC=10` - fp80 (extended precision). Native hardware support is available on x86 CPUs only. At present the only Fortran compiler that supports it is gfortran.
* `PREC=16` - fp128 (quadruple precision). This precision has no native hardware support on modern CPU architectures (except IBM Power) but it is supported by some compilers in the form of emulation.

Because calculations with fp128 floating-point numbers (even if compilers support them) use emulation on pretty much any modern CPU, they are very slow. As a rule of thumb, the speed of computations and memory allocation requirements of the codes in this collection on x86 CPUs scale as follows (the CPU time factors below are rough estimates; the exact scaling may vary depending on several aspects):

| Precision | CPU time (relative) | Memory amount (relative) |
| :---: | :---: | :---: |
| fp64 | x 1 | x 1 |
| fp80 | x 5 | x 2 |
| fp128 | x 100 | x 2 |

Therefore, quadruple precision (which is roughly a factor of ~100 slower than double precision) should be used only in special cases where it is absolutely necessary and feasible.

### Compilers

All ECGPACK codes are written in standard Fortran, use standard MPI calls and therefore, in principle, should be portable to any hardware and can be built with any Fortran compiler and MPI library. At present, however, they are routinely run and tested only on x86 hardware in a Linux environment. The following compilers/toolchains/precisions are currently supported and included in Makefiles:

| Toolchain kind | Compiler | MPI wrapper | MPI | Precision (real kind) | Allows MPI parallelism |
| :---: | :---: | :---: | :---: | :---: | :---: |
| GNU | gfortran | mpif90 | OpenMPI | 8  | yes |
| GNU | gfortran | mpif90 | OpenMPI | 10 | yes |
| GNU | gfortran | mpif90 | OpenMPI | 16 | no |
| Intel | ifort | mpiifort | Intel MPI | 8 | yes |
| Intel | ifort | mpiifort | Intel MPI | 16 | yes |
| Intel | ifx | mpiifx| Intel MPI | 8 | yes |
| Intel | ifx | mpiifx | Intel MPI | 16 | yes |
| Nvidia | nvfortran | mpif90 | CUDA-aware MPI (based on OpenMPI) | 8 | yes |

Note that at present the codes compiled with gfortran that use quadruple precision can only be executed in serial mode. This is a limitation of OpenMPI. A workaround is possible but has not been implemented yet.

To use Fortran compilers and MPI wrappers other than those listed above, the user needs to modify the Makefile accordingly.

### Linear algebra libraries

All energy codes (`CG_0S`, `RG_0S`, `RG_1P`, `RG_2D`, `RG_2P`) as an option can invoke the slow eigenvalue solver from standard LAPACK. This requires linking of these codes against some BLAS/LAPACK library. The choice of this library can be controlled with the `LINALG` argument of the `make` command:

* `LINALG=netlib` (default) - The bundled, lightly modified netlib reference BLAS/LAPACK subroutines are compiled from source (files `src/BLAS.f` and `src/LAPACK.f`). No external LAPACK/BLAS library is required in this case.
* `LINALG=mkl` - Intel Math Kernel Library (MKL).
* `LINALG=lblas` - An optimized BLAS/LAPACK exposed through the `-llapack -lblas` symbolic links. For example, if an NVHPC module is loaded last, this will lead to linking against BLAS/LAPACK that comes with NVIDIA HPC SDK. If an Easybuild's module `foss/2025b` is loaded last it will result in routing through FlexiBLAS, which will back-end directly into OpenBLAS.
* `LINALG=openblas` - OpenBLAS library is used through the `-lopenblas` link flag.
* `LINALG=aocl` - AMD Optimizing CPU Libraries (AOCL-BLAS and AOCL-LAPACK).

Note that an optimized BLAS/LAPACK library can be used only when `PREC=8`. For `PREC=10` and `PREC=16` only the `LINALG=netlib` option is available, because vendor-provided optimized BLAS/LAPACK libraries do not support extended or quadruple precision. In the off-diagonal matrix-element codes the `LINALG` argument is accepted for convenience in the `make` command but it has no effect because those codes do not use and do not link any BLAS/LAPACK library.

### Number of particles

It is very important to keep in mind that the number of particles is hardcoded at compile-time rather than passed as a runtime argument. An executable compiled for a specific number of particles will fail to run if the input file specifies a different count. This design constraint is strictly enforced for performance optimization, ensuring the code runs at maximum efficiency.

The number of particles is hardcoded in files `src/wp_def_*.f90` (here `*` stands for precision - 8, 10, or 16) by setting the value of parameter `Glob_AllowedNumOfParticles` accordingly. The batch-compile script `build.bash` does it automatically.

### Batch compilation of multiple code variants

The easiest way to compile all or some number of selected codes (basis type, number of particles, precisions, external libraries, etc.) in **one step** on a specific machine/OS using specific toolchains is to invoke the `build.bash` script located in the root directory. This script requires arguments. **Please read its source or run it with no arguments** to see instructions regarding how to run it properly. When this script is run, it will automatically move all individual binaries built to directory `ecgpack/bin/<toolchainname>/<configuration>`, where `<toolchainname>` (e.g. `systemdefault`) is the name of the toolchain specified and `<configuration>` can be either `debug` (slow and unoptimized binary suitable for debugging) or `release` (fast and optimized binary suitable for production work). The binary files will be named `<code_name>_N<particle_number>_P<precision>_<linalg>` (e.g. `RG_0S_N4_P8_netlib` - for `RG_0S` code, 4 particles, double precision, bundled netlib BLAS/LAPACK). The `<linalg>` suffix stands for the selected linear algebra library (`netlib`, `mkl`, `lblas`, `openblas`, or `aocl`) and is always present.

An example of executing the `build.bash` script to build production (optimized) binaries for `RG_0S` and `RG_1P` codes using double and extended precision and bundled BLAS/LAPACK source for the case of 4, 5, and 6 particles:

```bash
./build.bash machine=linux-generic toolchain=systemdefault config=release code=RG_0S,RG_1P nparticles=4,5,6 precision=8,10 linalg=netlib
```

Note that `systemdefault` toolchain assumes that the system's default `mpif90` wrapper is accessible out of the box without loading any environment modules - regardless of the underlying compiler or MPI implementation it wraps.

Script `build.bash` can use some common toolchains available in HPC systems/environments that are deployed with `Easybuild` - an open-source software management tool for scientific software, compilers, MPI libraries, BLAS/LAPACK libraries, and related packages. The following `Easybuild` toolchains can be invoked:

| Toolchain argument in `make` | Easybuild module | Description |
| :---: | :---: | :---: |
| systemdefault  | n/a | Default compiler and MPI on a given system referenced by `mpif90` |
| foss-2025b | foss/2025b | GNU Compiler Collection based compiler toolchain |
| foss-2025a | foss/2025a | GNU Compiler Collection based compiler toolchain |
| foss-2024a | foss/2024a | GNU Compiler Collection based compiler toolchain |
| foss-2023b | foss/2023b | GNU Compiler Collection based compiler toolchain |
| intel-2025b | intel/2025b | Intel compilers and libraries |
| intel-2025a | intel/2025a | Intel compilers and libraries |
| intel-2024a | intel/2024a | Intel compilers and libraries |  
| intel-2023b | intel/2023b | Intel compilers and libraries |
| nvhpc-25.9 | NVHPC/25.9-CUDA-12.9.1 | Nvidia compilers and libraries included in NVHPC SDK |
| nvhpc-25.3 | NVHPC/25.3-CUDA-12.8.0 | Nvidia compilers and libraries included in NVHPC SDK |

Users can certainly easily modify `build.bash` to suit their own environment.

## Execution

Any compiled binary file can be executed using multiple MPI processes in a standard way:

```bash
mpirun -np <NPROCS> <BINARYFILE>
```

For example, suppose a user has compiled all codes with a help of `build.bash` script (invoking the `systemdefault` toolchain, which on most Linux machines means `gfortran` compiler with `OpenMPI`). Then the user creates a work directory `ecgpack/jobs/test` for a test calculation of a 4-particle system with the `RG_0S` basis (states with $L=0$ angular momentum quantum number) and places a relevant input file `inout.txt` in that directory. To launch a calculation with double precision using 32 CPU cores the execution command will look as follows:

```bash
mpirun -np 32 ../../bin/systemdefault/release/RG_0S_N4_P8_netlib
```

Note that the input file `inout.txt` (as well as other relevant files, if any) should be located in the same directory where the execution of the binary file takes place.

## Note on parallel scaling

In ECG calculations, two major computational bottlenecks typically consume the most time:

1. Evaluation of matrix elements
2. Solving the generalized eigenvalue problem (typically by updating the solution obtained in the previous step of the calculation)

The distribution of time between these two tasks varies significantly depending on the calculation type, basis set, number of particles (especially identical ones), basis size, and network interconnect speed.

Which task dominates largely determines the parallel scalability of the simulation:

* Task #1 (Matrix Elements): Highly parallelizable ("embarrassingly parallel"), showing excellent scaling with the number of MPI processes.
* Task #2 (Eigenvalue Problem): Scales significantly worse due to communication overhead and other factors.

Consequently, parallel efficiency is highly system-dependent:

* Poor Scaling: Systems with a small number of particles (e.g., 3) and, at the same time, very large bases (thousands of functions). Here, Task #2 either dominates or takes a significant fraction of CPU time, leading to poor MPI scaling. In this case it makes little sense to use more than a few CPU cores in calculations.
* Linear Scaling: Systems with 5, 6, or more electrons (or other identical particles) and basis sizes that are not gigantic. Here, Task #1 dominates heavily, resulting in nearly linear scaling even across hundreds of MPI processes and multiple nodes, provided a fast network interconnect is used (e.g. Infiniband).

## Note on scaling with the number of particles

The single most important factor determining the cost of a calculation is the number of particles $N$, which is fixed at compile time (see [Number of particles](#number-of-particles) above). The table below gives rough relative timings for the evaluation of one complete set of matrix elements between a pair of basis functions, that is the Hamiltonian matrix element, the overlap matrix element, and their gradients with respect to the nonlinear parameters, as these are evaluated in wave function optimization runs. Permutational symmetry is *not* taken into account here (see the remarks below the table). All values are normalized so that the cheapest case, `RG_0S` with 3 particles, is equal to 1.

| $N$ | RG_0S | RG_1P | RG_2D | RG_2P |
| :---: | :---: | :---: | :---: | :---: |
| 3 | 1.0 | 1.4 | 2.1 | 2.0 |
| 4 | 2.3 | 3.0 | 4.7 | 4.8 |
| 5 | 4.8 | 5.9 | 7.9 | 7.9 |
| 6 | 8.6 | 11 | 14 | 14 |
| 7 | 14 | 16 | 20 | 19 |
| 8 | 21 | 24 | 29 | 28 |
| 9 | 21 | 25 | 33 | 33 |

A few remarks on these numbers and how they should be read:

* The entries are a weighted mix of the two kinds of calls made by the optimizer. During a wave function optimization the minimizer requests the energy alone roughly twice as often as it requests the energy together with the gradient, so the numbers above correspond to a 2:1 ratio of such calls. Gradient evaluations are nevertheless what dominates the total time, because a call with gradients is several times more expensive than a call without them.
* The timings were measured in double precision (fp64) on a single core of one x86 node with the GNU toolchain. Timings obtained on other machines and with other compilers may deviate from these by tens of percent, but the relative pattern, both across basis types and across $N$, is quite stable.
* Bases with a higher angular momentum are more expensive: at a fixed $N$, `RG_1P` costs about 15–40% more than `RG_0S`, while `RG_2D` and `RG_2P` cost anywhere from about 40% more to about twice as much, with the largest relative penalty at small $N$. `RG_2D` and `RG_2P` are practically indistinguishable from each other in this respect, and the small differences between their columns above are well within the noise of the measurement.
* Permutational symmetry acts as an additional multiplicative factor on top of the numbers in the table. The complete matrix element is a sum over the terms retained in the symmetry projector ( $Y^\dagger Y$ operator), and each such term costs about as much as one entry of the table. The number of these terms grows very rapidly, normally as $\sim n!$ , where $n$ is the number of identical particles. For systems with many identical particles it easily exceeds the cost increase due to $N$ itself. This is the reason why calculations with 5 or more electrons are heavily dominated by the evaluation of matrix elements (see [Note on parallel scaling](#note-on-parallel-scaling) above).
* Finally, the total number of matrix elements to be evaluated grows as $K(K+1)/2$ with the basis size $K$, so the total time spent per optimization step is the product of that number, the cost from the table, and the permutational factor.

