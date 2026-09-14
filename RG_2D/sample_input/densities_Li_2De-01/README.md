# Lithium ${}^2D^e$ density example

This directory contains a 100-function `RG_2D` example for the four-particle Lithium atom (Li-7 isotope). It evaluates coordinate- and momentum-space particle densities and pair correlation functions on a piecewise cylindrical grid.

## Build and run

To reproduce the results, from the repository root, build the double-precision, four-particle executable:

```bash
./build.bash machine=linux-generic toolchain=systemdefault config=release code=RG_2D nparticles=4 precision=8 linalg=netlib
```

then generate the grid and run the ECGPACK calculations:

```bash
cd RG_2D/sample_input/densities_Li_2De-01
python3 create_cyl_grid.py
mkdir -p output_cyl
mpirun -np 1 ../../../bin/systemdefault/release/RG_2D_N4_P8_netlib
mv -f cf.dat dens.dat mcf.dat mdens.dat expvals.txt spinData.txt swapfile.dat output_cyl/
```

The paths in this command sequence are relative to the working directory `RG_2D/sample_input/densities_Li_2De-01`, and the example starts one MPI process. To change the number of MPI processes, edit:

```bash
mpirun -np 1 ../../../bin/systemdefault/release/RG_2D_N4_P8_netlib
```

The final command moves the generated files into `output_cyl/`, overwriting the sample outputs. Remove that command if you want to keep the files next to `inout.txt`.

On a four-core Intel Core i5-1135G7 system under WSL, the supplied 16,836-point grid and 100-function basis take approximately 30 seconds with four MPI processes. This is a separate four-process benchmark; the command above uses one process. Build time is not included. Runtime varies with the compiler, linear-algebra library, processor, system load, basis size, grid size, and MPI process count.

The executable must be built for the `PARTICLES 4` value in `inout.txt`.

## Density commands in `inout.txt`

The general command forms are:

```text
DENSITIES <G|I|Q> <basis_size> <cf_grid> <cf_output> <dens_grid> <dens_output>
MOMT_DENS <G|I|Q> <basis_size> <mom_cf_grid> <mom_cf_output> <mom_dens_grid> <mom_dens_output>
```

`G`, `I`, and `Q` select the direct generalized solver, inverse iteration, and QR-based inverse iteration, respectively; this example uses inverse iteration (`I`). `<basis_size>` must equal the current basis size. The first file pair controls correlation functions and the second controls particle densities. Use `none` as the corresponding grid name to disable one family, while retaining all positional arguments. The same grid may be supplied to both families, as in this example.

`DENSITIES` evaluates center-of-mass-frame coordinate-space functions. `MOMT_DENS` evaluates their momentum-space analogues. Both commands also perform the `EXPC_VALS` work, so no separate `EXPC_VALS` command is required.

The BBOP commands create the following tables in the working directory; the run instructions above then archive them in `output_cyl/`:

- `output_cyl/cf.dat`: coordinate-space correlation functions;
- `output_cyl/dens.dat`: coordinate-space particle densities;
- `output_cyl/mcf.dat`: momentum-space correlation functions;
- `output_cyl/mdens.dat`: momentum-space particle densities.

Each table starts with a `#` header. Coordinate-space tables use cylindrical variables $(\xi_\rho,\xi_z)$, correlation functions $g$, and densities $\rho$. Momentum-space tables use $(\eta_\rho,\eta_z)$, correlation functions $f$, and densities $\varrho$. `expvals.txt`, `spinData.txt`, and `swapfile.dat` are ancillary calculation files.

## Grid generation

`create_cyl_grid.py` constructs `sample_cyl_ecg_grid.dat` as the tensor product of separately generated $\rho$ and $z$ axes. Edit `axis_values`, `axis_powers`, `rho_counts`, and `z_counts` in `main()` to change the piecewise grid. Counts are numbers of intervals, shared boundaries are emitted once, $\rho$ is nonnegative, and the positive $z$ axis is reflected to form an increasing signed axis. The supplied grid contains 92 $\rho$ points and 183 $z$ points, for 16,836 points in total, and extends to 48 a.u.

The same two columns are interpreted as $(\xi_\rho,\xi_z)$ by `DENSITIES` and as $(\eta_\rho,\eta_z)$ by `MOMT_DENS`; all quantities are in atomic units. Re-running the script replaces `sample_cyl_ecg_grid.dat`, and re-running ECGPACK replaces the four requested tables. `create_radial_grid.py` remains available as an alternative generator, but no radial-grid output tables or plots are supplied in this example.

## Plots

The four gnuplot scripts read the cylindrical-grid tables in `output_cyl/` and write direct functions and axisymmetric, volume-weighted distributions to `cyl_plots/`:

```bash
gnuplot create_plots_coor_cf.gp
gnuplot create_plots_coor_dens.gp
gnuplot create_plots_mom_cf.gp
gnuplot create_plots_mom_dens.gp
```

Only cylindrical-grid sample plots are provided. To read tables from another directory or write plots elsewhere, define the `output_dir` or `plot_dir` gnuplot variables.

The plotting windows intentionally emphasize the significant, rapidly varying part of each function rather than its long, low-valued tail. The committed examples are in `cyl_plots/`.

## Nuclear-density resolution

Particle 1 is the finite-mass lithium nucleus. Its center-of-mass-frame density is extremely narrow compared with the electronic density, with a three-orders-of-magnitude difference in width. Thus, a highly resolved region within approximately $10^{-3}$ a.u. of the origin is necessary if the nuclear density is of interest.

## Pointwise density and cylindrical distribution

The direct density and the plotted cylindrical distribution are different quantities. The axisymmetric volume element contributes $2\pi\xi_{\rho}$, so $2\pi\xi_{\rho}\,\rho_i(\xi_{\rho},\xi_z)$ vanishes on the symmetry axis even when the pointwise density does not.
