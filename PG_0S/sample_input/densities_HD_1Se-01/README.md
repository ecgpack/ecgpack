# HD ${}^1S^e$ correlation function example

This directory contains a 10-function `PG_0S` example for the four-particle HD molecule (deuteron, proton and two electrons, all treated non-adiabatically). It evaluates the nucleus-nucleus pair correlation function on the supplied radial grid. In `PG_0S` this is the only correlation function available; particle densities and the momentum-space functions of `MOMT_DENS` are not available.

## Build and run

To reproduce the results, from the repository root, build the double-precision, four-particle executable:

```bash
./build.bash machine=linux-generic toolchain=systemdefault config=release code=PG_0S nparticles=4 precision=8 linalg=netlib
```

then generate the grid and run the ECGPACK calculation:

```bash
cd PG_0S/sample_input/densities_HD_1Se-01
python3 create_grid.py
mkdir -p output_cyl
mpirun -np 1 ../../../bin/systemdefault/release/PG_0S_N4_P8_netlib
mv -f cf.dat expvals.txt swapfile.dat output_cyl/
```

The paths in this command sequence are relative to the working directory `PG_0S/sample_input/densities_HD_1Se-01`, and the example starts one MPI process. To change the number of MPI processes, edit:

```bash
mpirun -np 1 ../../../bin/systemdefault/release/PG_0S_N4_P8_netlib
```

The final command moves the generated files into `output_cyl/`. Remove that command if you want to keep the files next to `inout.txt`. `PG_0S` uses a one-dimensional radial grid; the directory name `output_cyl` is kept for consistency with the density examples of the other codes.

On an Intel Core i9-7900X system under WSL, the supplied 601-point grid and 10-function basis take less than 1 second with one MPI process (release build, gfortran). Build time is not included.

The executable must be built for the `PARTICLES 4` value in `inout.txt`.

## Density commands in `inout.txt`

The general command form is the same as in the other codes:

```text
DENSITIES <G|I|Q> <basis_size> <cf_grid> <cf_output> <dens_grid> <dens_output>
```

`G`, `I`, and `Q` select the direct generalized solver, inverse iteration, and QR-based inverse iteration, respectively; this example uses the direct solver (`G`). `<basis_size>` must equal the current basis size. The first file pair controls the correlation function and the second controls particle densities. Particle densities are not available in `PG_0S`, so keep `none none` for the second pair; a density grid other than `none` is ignored with warning WC0195. Use `none` as the correlation-function grid name to skip the correlation function. The command used in this example is:

```text
DENSITIES G 10 cf_grid.dat cf.dat none none
```

`DENSITIES` evaluates the center-of-mass-frame coordinate-space correlation function. It also performs the `EXPC_VALS` work, so no separate `EXPC_VALS` command is required.

The command creates the following table in the working directory; the run instructions above then archive it in `output_cyl/`:

- `output_cyl/cf.dat`: the nucleus-nucleus correlation function.

The table starts with the `#` header `#r g1`, followed by one line per grid point with the radius $r$ and $g_1(r)$. $g_1$ refers to particles 0 and 1, the two nuclei, whose distance is $r_1$. `expvals.txt` and `swapfile.dat` are ancillary calculation files.

## Grid generation

`create_grid.py` constructs the one-column radial file `cf_grid.dat`. The variables `counts`, `powers`, and `x_values` in `main()` can be edited to change the number and placement of points. Counts are numbers of intervals, shared segment boundaries are emitted once, and `powers` control within-segment clustering. The supplied grid is one segment of 600 intervals from 0 to 6 bohr (601 points, step 0.01).

All quantities are in atomic units. Re-running the script replaces `cf_grid.dat`, and re-running ECGPACK replaces `cf.dat`.

## Plots

No plotting scripts are supplied. The correlation function and the radial distribution can be viewed with gnuplot:

```bash
gnuplot -p -e 'set y2tics; plot "output_cyl/cf.dat" using 1:2 axes x1y2 with lines title "g1 (right axis)", "" using 1:(4*pi*$1**2*$2) with lines title "4 pi r^2 g1"'
```

## Grid resolution

For this basis, $g_1$ is concentrated around the equilibrium internuclear distance: it falls to half of its maximum at 1.05 and 1.68 bohr, and 99.9% of the probability lies between 0.66 and 2.69 bohr. The uniform step of 0.01 bohr resolves it well. For other bases, place the points where $g_1$ is large rather than extending the grid.

As a check, $\int 4\pi r^2 g_1(r)\,dr$ over the grid should be 1, and $\int 4\pi r^3 g_1\,dr$ and $\int 4\pi r^4 g_1\,dr$ should reproduce $\langle r_1\rangle=1.4787864054$ and $\langle r_1^2\rangle=2.2625366626$ printed by the same run. With the supplied grid they agree to better than $10^{-14}$.

## Pointwise correlation function and radial distribution

The pointwise function $g_1(r)$ and the radial distribution $4\pi r^2 g_1(r)$ are different quantities. The radial distribution is the probability density of the internuclear distance; for this basis it peaks at 1.47 bohr, while $g_1$ itself peaks at 1.37 bohr with the value 0.05992. Because every basis function contains the factor $r_1^{2 m_k}$ with $m_k \ge 1$, $g_1(0)=0$: the two nuclei never coincide.
