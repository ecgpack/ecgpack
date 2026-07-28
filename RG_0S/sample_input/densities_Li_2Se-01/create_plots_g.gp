#!/usr/bin/env gnuplot
#
# Gnuplot script: plot correlation functions from cf.dat
# Column 1: r, columns 2..N: g-functions named in the header line,
# e.g.  g1, g2, g3, g12, g13, g23  (g<i> one-particle, g<ij> pair).
# Works for an arbitrary number / naming of g-columns.
#
# For each column two plots are produced: g(r) itself and the radial
# distribution 4 pi r^2 g(r). For each curve the upper end of the
# x-range is trimmed automatically: the plot extends only up to the
# largest r at which the curve is still above TAIL_FRAC of its own
# maximum value (i.e. the decayed tail is cut).

set terminal pngcairo enhanced size 2000,1400 font "Helvetica,28"

set xlabel "r"
set grid

TAIL_FRAC = 0.005   # 0.5% of the peak value

# Column names taken from the header line (fields 2..N, leading '#' stripped)
names = system("awk 'NR==1{for(i=2;i<=NF;i++){n=$i; sub(/^#/,\"\",n); printf \"%s \", n}}' cf.dat")

# Grid lower bound
stats "cf.dat" using 1 nooutput
rmin = STATS_min

ng = words(names)

do for [i=1:ng] {
    col  = i + 1
    name = word(names, i)
    sub  = name[2:strlen(name)]   # subscript = characters after the leading 'g'

    # ---- g(r) ----
    stats "cf.dat" using col nooutput
    thr = TAIL_FRAC * STATS_max
    stats "cf.dat" using (column(col) >= thr ? $1 : 1/0) nooutput
    set xrange [rmin:STATS_max]

    set output sprintf("%s.png", name)
    set ylabel sprintf("g_{%s}(r)", sub)
    set title  sprintf("Correlation function g_{%s}(r)", sub)
    plot "cf.dat" using 1:col with lines lw 4 \
         title sprintf("g_{%s}(r)", sub)

    # ---- 4 pi r^2 g(r) ----
    stats "cf.dat" using (4*pi*$1**2*column(col)) nooutput
    thr = TAIL_FRAC * STATS_max
    stats "cf.dat" using ((4*pi*$1**2*column(col)) >= thr ? $1 : 1/0) nooutput
    set xrange [rmin:STATS_max]

    set output sprintf("%s_4pir2.png", name)
    set ylabel sprintf("4{/Symbol p} r^2 g_{%s}(r)", sub)
    set title  sprintf("Radial distribution 4{/Symbol p} r^2 g_{%s}(r)", sub)
    plot "cf.dat" using 1:(4*pi*$1**2*column(col)) with lines lw 4 \
         title sprintf("4{/Symbol p} r^2 g_{%s}(r)", sub)
}
