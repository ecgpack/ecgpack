# Lithium $2\,{}^2S^e$ density example

This directory contains a 100-function `RG_0S` example for the four-particle lithium atom. It evaluates coordinate- and momentum-space particle densities and pair correlation functions on the supplied radial grid.

## Build and run

To reproduce the results, from the repository root, build the double-precision, four-particle executable:

```bash
./build.bash machine=linux-generic toolchain=systemdefault config=release code=RG_0S nparticles=4 precision=8 linalg=netlib
```

then generate the grid and run the ECGPACK calculations:

```bash
cd RG_0S/sample_input/densities_Li_2Se-01
python3 create_grid.py
mkdir -p output_cyl
mpirun -np 1 ../../../bin/systemdefault/release/RG_0S_N4_P8_netlib
mv -f cf.dat dens.dat mcf.dat mdens.dat expvals.txt spinData.txt swapfile.dat output_cyl/
```

The paths in this command sequence are relative to the working directory `RG_0S/sample_input/densities_Li_2Se-01`, and the example starts one MPI process. To change the number of MPI processes, edit:

```bash
mpirun -np 1 ../../../bin/systemdefault/release/RG_0S_N4_P8_netlib
```

The final command moves the generated files into `output_cyl/`, overwriting the sample outputs. Remove that command if you want to keep the files next to `inout.txt`. RG_0S uses a one-dimensional radial grid; the common `output_cyl` directory name is retained across all four density examples.

On a four-core Intel Core i5-1135G7 system under WSL, the supplied 1,251-point grid and 100-function basis take approximately 3 seconds with four MPI processes. This is a separate four-process benchmark; the command above uses one process. Build time is not included. Runtime varies with the compiler, linear-algebra library, processor, system load, basis size, grid size, and MPI process count.

The executable must be built for the `PARTICLES 4` value in `inout.txt`.

## Density commands in `inout.txt`

The general command forms are:

```text
DENSITIES <G|I> <basis_size> <cf_grid> <cf_output> <dens_grid> <dens_output>
MOMT_DENS <G|I> <basis_size> <mom_cf_grid> <mom_cf_output> <mom_dens_grid> <mom_dens_output>
```

`G` and `I` select the generalized-eigenvalue solver; this example uses inverse iteration (`I`). `<basis_size>` must equal the current basis size. The first file pair controls correlation functions and the second controls particle densities. Use `none` as the corresponding grid name to disable one family, while retaining all positional arguments. The same grid may be supplied to both families, as in this example.

`DENSITIES` evaluates center-of-mass-frame coordinate-space functions. `MOMT_DENS` evaluates their momentum-space analogues. Both commands also perform the `EXPC_VALS` work, so no separate `EXPC_VALS` command is required.

The BBOP commands create the following tables in the working directory; the run instructions above then archive them in `output_cyl/`:

- `output_cyl/cf.dat`: coordinate-space correlation functions;
- `output_cyl/dens.dat`: coordinate-space particle densities;
- `output_cyl/mcf.dat`: momentum-space correlation functions;
- `output_cyl/mdens.dat`: momentum-space particle densities.

Each table starts with a `#` header. Coordinate-space tables use radius $r$, correlation functions $g$, and densities $\rho$. Momentum-space tables use momentum variable $\eta$, correlation functions $f$, and densities $\varrho$. `expvals.txt`, `spinData.txt`, and `swapfile.dat` are ancillary calculation files.

## Grid generation

`create_grid.py` constructs the one-column radial file `sample_radial_ecg_grid.dat`. The variables `counts`, `powers`, and `x_values` in `main()` can be edited to change the number and placement of points. Counts are numbers of intervals, shared segment boundaries are emitted once, and `powers` control within-segment clustering. The initial 0 to 0.002 segment uses 200 intervals (201 endpoint-inclusive samples) to resolve the nuclear density, and the grid ends at 64 a.u.

The same numerical grid is interpreted as a coordinate radius by `DENSITIES` and as a momentum magnitude by `MOMT_DENS`; all quantities are in atomic units. Re-running the script replaces `sample_radial_ecg_grid.dat`, and re-running ECGPACK replaces the four requested tables.

## Plots

The four gnuplot scripts read the sample tables in `output_cyl/` and write direct functions and volume-weighted radial distributions to `cyl_plots/`:

```bash
gnuplot create_plots_coor_cf.gp
gnuplot create_plots_coor_dens.gp
gnuplot create_plots_mom_cf.gp
gnuplot create_plots_mom_dens.gp
```

Only this sample result set is provided. To read tables from another directory or write plots elsewhere, define the `output_dir` or `plot_dir` gnuplot variables.

The plotting windows intentionally emphasize the significant, rapidly varying part of each function rather than its long, low-valued tail. The committed examples are in `cyl_plots/`.

## Nuclear-density resolution

Particle 1 is the finite-mass lithium nucleus. Its center-of-mass-frame density is extremely narrow compared with the electronic density, with a three-orders-of-magnitude difference in width. Thus, a highly resolved region within approximately $10^{-3}$ a.u. of the origin is necessary if the nuclear density is of interest.

## Pointwise density and radial distribution

The pointwise density and the plotted radial distribution are different quantities. The spherical volume element contributes $4\pi r^2$, so $4\pi r^2\rho_i(r)$ vanishes at the origin even when the pointwise density $\rho_i(0)$ does not.
