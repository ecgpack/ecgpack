#!/usr/bin/env gnuplot
#
# Plot radial correlation functions from cf.dat.

set terminal pngcairo enhanced size 2000,1400 font "Helvetica,28"
set xlabel "r"
set grid
TAIL_FRAC = 0.01
if (!exists("output_dir")) output_dir = "output_cyl"
if (!exists("plot_dir")) plot_dir = "cyl_plots"
system(sprintf("mkdir -p %s", plot_dir))
datafile = output_dir . "/cf.dat"
cleanfile = "cf_plot.tmp"

# Normalize exponent fields that omit the Fortran E marker.
system(sprintf("awk '{for(i=1;i<=NF;i++) if ($i !~ /[EeDd]/ && $i ~ /[0-9][+-][0-9][0-9][0-9]$/) $i=substr($i,1,length($i)-4) \"E\" substr($i,length($i)-3); print}' %s > %s",datafile,cleanfile))

# Read the function names directly from the output-table header.
names = system(sprintf("awk 'NR==1{for(i=2;i<=NF;i++){n=$i; sub(/^#/,\"\",n); printf \"%%s \",n}}' %s", datafile))
stats cleanfile using 1 nooutput
rmin = STATS_min

# Plot each function and its radial integration kernel.
do for [i=1:words(names)] {
    col = i + 1
    name = word(names,i)
    subscript = name[2:strlen(name)]
    stats cleanfile using col nooutput
    thr = TAIL_FRAC * STATS_max
    stats cleanfile using (column(col) >= thr ? $1 : 1/0) nooutput
    set xrange [rmin:STATS_max]
    set output sprintf("%s/%s.png", plot_dir, name)
    set ylabel sprintf("g_{%s}(r)",subscript)
    set title sprintf("Correlation function g_{%s}(r)",subscript)
    plot cleanfile using 1:col with lines lw 4 title sprintf("g_{%s}(r)",subscript)

    stats cleanfile using (4*pi*$1**2*column(col)) nooutput
    thr = TAIL_FRAC * STATS_max
    stats cleanfile using ((4*pi*$1**2*column(col)) >= thr ? $1 : 1/0) nooutput
    set xrange [rmin:STATS_max]
    set output sprintf("%s/%s_4pir2.png", plot_dir, name)
    set ylabel sprintf("4{/Symbol p}r^2 g_{%s}(r)",subscript)
    set title sprintf("Radial distribution 4{/Symbol p}r^2 g_{%s}(r)",subscript)
    plot cleanfile using 1:(4*pi*$1**2*column(col)) with lines lw 4 title sprintf("4{/Symbol p}r^2 g_{%s}(r)",subscript)
}
system(sprintf("rm -f %s",cleanfile))
