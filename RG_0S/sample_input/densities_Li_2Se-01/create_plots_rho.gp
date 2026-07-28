#!/usr/bin/env gnuplot
#
# Gnuplot script: plot radial densities from dens.dat
# Column 1: r, columns 2..N: rho1(r), rho2(r), ...
# Works for an arbitrary number of density columns.
#
# For each curve the upper end of the x-range is trimmed automatically:
# the plot extends only up to the largest r at which the curve is still
# above TAIL_FRAC of its own maximum value (i.e. the decayed tail is cut).

set terminal pngcairo enhanced size 2000,1400 font "Helvetica,28"

set xlabel "r"
set grid

TAIL_FRAC = 0.005   # 0.5% of the peak value

# Determine the grid range and the number of columns in the data file
stats "dens.dat" using 1 nooutput
rmin  = STATS_min
ndens = STATS_columns - 1   # number of density columns

do for [i=1:ndens] {
    col = i + 1

    # ---- rho_i(r) ----
    stats "dens.dat" using col nooutput
    thr = TAIL_FRAC * STATS_max
    stats "dens.dat" using (column(col) >= thr ? $1 : 1/0) nooutput
    set xrange [rmin:STATS_max]

    set output sprintf("rho%d.png", i)
    set ylabel sprintf("{/Symbol r}_%d(r)", i)
    set title  sprintf("Density {/Symbol r}_%d(r)", i)
    plot "dens.dat" using 1:col with lines lw 4 \
         title sprintf("{/Symbol r}_%d(r)", i)

    # ---- 4 pi r^2 rho_i(r) ----
    stats "dens.dat" using (4*pi*$1**2*column(col)) nooutput
    thr = TAIL_FRAC * STATS_max
    stats "dens.dat" using ((4*pi*$1**2*column(col)) >= thr ? $1 : 1/0) nooutput
    set xrange [rmin:STATS_max]

    set output sprintf("rho%d_4pir2.png", i)
    set ylabel sprintf("4{/Symbol p} r^2 {/Symbol r}_%d(r)", i)
    set title  sprintf("Radial distribution 4{/Symbol p} r^2 {/Symbol r}_%d(r)", i)
    plot "dens.dat" using 1:(4*pi*$1**2*column(col)) with lines lw 4 \
         title sprintf("4{/Symbol p} r^2 {/Symbol r}_%d(r)", i)
}
