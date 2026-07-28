# Description

Sample file `inout.txt` in this directory contains a basis of 100 ECG functions for the lowest $^2 P^o$ state of Li atom (Li-7 isotope) and sets up calculations of particle distributions in the coordinate space: particle densities in the center-of-mass frame, $\rho_i(\xi_\rho,\xi_z)$ , and pair correlation functions $g_{i}(\xi_\rho,\xi_z)$ and $g_{ij}(\xi_\rho,\xi_z)$ . Here $\xi_\rho=\sqrt{\xi_x^2+\xi_y^2}$ . Expectation values of various operators are also computed in one go.
The calculations with this input file `inout.txt` and the grids used in the example should take roughly 15 seconds on a single CPU core when double precision is used (gfortran compiler / AMD Ryzen AI 9 HX 370).

The instruction that is used in `inout.txt` to run calculations of distributions has the following format:  
`DENSITIES I 100 cf_grid.dat cf.dat dens_grid.dat dens.dat`  
Here `cf_grid.dat` and `dens_grid.dat` are files that contain the 2D grids of $\xi_\rho$ and $\xi_z$ values where the correlation functions and densities need to be computed. Both of these grid files have two columns and look like this:  

```
0.000000 -1.500000
0.000000 -1.480000
0.000000 -1.460000
0.000000 -1.440000
...
...
...
```

In principle, the grids for both the correlation functions and densities could come from the same file (e.g. `grid.dat`) if the expected range of all distributions is comparable. The grids do not need to be uniform. The points do not have to be sorted in any particular order, although doing so may be convenient for further processing and visualization later.

For atomic systems with a finite-mass nucleus, such as Li, the density of the first particle (nucleus) is expected to be highly localized around the origin (center of mass), $\xi=0$ . In other words, the effective range of the nuclear distribution is about 3 orders of magnitude shorter than that for electrons. If both the electronic and nuclear densities are needed then the calculations for them should be run independently and use different grids.  

When the calculations are executed the grid files must be in the same directory where the `inout.txt` is located and where the program runs.

A Bash script `create_grid.bash` is also provided for convenience. One can adjust grid parameters in this script. When executed it creates files `cf_grid.dat` and `dens_grid.dat`.

After correlation functions and densities have been calculated, files `cf.dat` and `dens.dat` are created in the execution directory. They have the following format:

```
      #xi_rho                     xi_z                     g1                      g12
  0.0000000000000000E+00 -0.1500000000000000E+01  0.3963529627411734E-02  0.6918869335184548E-02
  0.0000000000000000E+00 -0.1480000000000000E+01  0.4132513711459111E-02  0.7148840535921017E-02
  0.0000000000000000E+00 -0.1460000000000000E+01  0.4318856778477886E-02  0.7397359829233206E-02
  0.0000000000000000E+00 -0.1440000000000000E+01  0.4524441074520022E-02  0.7665890303721264E-02
  ...
  ...
  ...
```

```
      #xi_rho                     xi_z                    rho1                    rho2
  0.0000000000000000E+00 -0.1500000000000000E+01  0.0000000000000000E+00  0.3963670099161181E-02
  0.0000000000000000E+00 -0.1480000000000000E+01  0.0000000000000000E+00  0.4132621704230993E-02
  0.0000000000000000E+00 -0.1460000000000000E+01  0.0000000000000000E+00  0.4318929489442429E-02
  0.0000000000000000E+00 -0.1440000000000000E+01  0.0000000000000000E+00  0.4524475524150619E-02
  ...
  ...
  ...
```

Note that if there are identical particles in the system (e.g. three electrons in Li atom), the number of $\rho_{i}$ distributions written in the `dens.dat` file is less than the number of particles. In the case of Li atom $\rho_1(\xi_\rho,\xi_z)$ is the nuclear density, while $\rho_2(\xi_\rho,\xi_z)=\rho_3(\xi_\rho,\xi_z)=\rho_4(\xi_\rho,\xi_z)$ are the electronic densities (so only $\rho_2(\xi_\rho,\xi_z)$ is written). Likewise, the number of pair correlation functions $g_i(\xi_\rho,\xi_z)$ and $g_{ij}(\xi_\rho,\xi_z)$ is less than the number of all possible pairs. In the case of Li, the only physically different pair correlation functions are the nucleus-electron one, $g_1=g_2=g_3$, and the electron-electron one, $g_{12}=g_{13}=g_{23}$.

The data in `cf.dat` and `dens.dat` can be further processed or visualized according to the needs of the user.

As an example of further processing, two gnuplot scripts (`create_plots_g.gp` and `create_plots_rho.gp`) are provided. They can be executed as follows:

```bash
gnuplot create_plots_g.gp
gnuplot create_plots_rho.gp
```

They read files `cf.dat` and `dens.dat`, respectively, and create contour plots in the `png` format for all distributions. Two plots are created for each 2D distribution: the function itself, e.g. $g_{12}(\xi_\rho,\xi_z)$ or $\rho_1(\xi_\rho,\xi_z)$ , and the corresponding axial probability density, e.g. $2 \pi \xi_\rho g_{12}(\xi_\rho,\xi_z)$ or $2 \pi \xi_\rho \rho_1(\xi_\rho,\xi_z)$ . The resulting `png` files generated by gnuplot for the considered case of Li atom are provided.  

![electron density in the center of mass frame](rho2.png)  

![electron radial probability in the center of mass frame](rho2_axial_dist.png)  

![nucleus-electron correlation function](g1.png)  

![electron-electron correlation function](g12.png)  

![nucleus-electron radial probability distribution](g1_axial_dist.png)  

![electron-electron radial probability distribution](g12_axial_dist.png)