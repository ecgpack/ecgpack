#!/usr/bin/env gnuplot
#
# Plot radial particle densities from dens.dat.  This is the original
# sample script with a guard for distributions narrower than the sample grid.

set terminal pngcairo enhanced size 2000,1400 font "Helvetica,28"
set xlabel "r"
set grid
TAIL_FRAC = 0.01
if (!exists("output_dir")) output_dir = "output_cyl"
if (!exists("plot_dir")) plot_dir = "cyl_plots"
system(sprintf("mkdir -p %s", plot_dir))
datafile = output_dir . "/dens.dat"
cleanfile = "dens_plot.tmp"

# Normalize exponent fields that omit the Fortran E marker.
system(sprintf("awk '{for(i=1;i<=NF;i++) if ($i !~ /[EeDd]/ && $i ~ /[0-9][+-][0-9][0-9][0-9]$/) $i=substr($i,1,length($i)-4) \"E\" substr($i,length($i)-3); print}' %s > %s",datafile,cleanfile))

# Discover the grid extent and number of density columns.
stats cleanfile using 1 nooutput
rmin = STATS_min
ndens = STATS_columns - 1

# Plot each density and its radial integration kernel.
do for [i=1:ndens] {
    col = i + 1
    set autoscale y
    stats cleanfile using col nooutput
    thr = TAIL_FRAC * STATS_max
    stats cleanfile using (column(col) >= thr ? $1 : 1/0) nooutput
    xmax = STATS_max
    if (xmax <= rmin) { xmax = 0.1 }
    set xrange [rmin:xmax]
    set output sprintf("%s/rho%d.png", plot_dir, i)
    set ylabel sprintf("{/Symbol r}_{%d}(r)", i)
    set title sprintf("Particle density {/Symbol r}_{%d}(r)", i)
    plot cleanfile using 1:col with lines lw 4 title sprintf("{/Symbol r}_{%d}(r)", i)

    set autoscale y
    stats cleanfile using (4*pi*$1**2*column(col)) nooutput
    peak = STATS_max
    if (peak > 0) {
        thr = TAIL_FRAC * peak
        stats cleanfile using ((4*pi*$1**2*column(col)) >= thr ? $1 : 1/0) nooutput
        xmax = STATS_max
    } else {
        xmax = 0.1
        set yrange [0:1]
    }
    set xrange [rmin:xmax]
    set output sprintf("%s/rho%d_4pir2.png", plot_dir, i)
    set ylabel sprintf("4{/Symbol p}r^2 {/Symbol r}_{%d}(r)", i)
    set title sprintf("Radial distribution 4{/Symbol p}r^2 {/Symbol r}_{%d}(r)", i)
    plot cleanfile using 1:(4*pi*$1**2*column(col)) with lines lw 4 title sprintf("4{/Symbol p}r^2 {/Symbol r}_{%d}(r)", i)
}
system(sprintf("rm -f %s",cleanfile))
