#!/usr/bin/env gnuplot
# ---------------------------------------------------------------------------
# create_plots_coor_dens.gp
#
# Reads dens.dat:
#   col 1 = xi_rho   (x-axis)
#   col 2 = xi_z     (y-axis)
#   col 3..N = densities rho_i (names taken from the header: rho1, rho2, ...)
#
# For every data column it produces two PNG contour/density plots:
#   <name>.png          ->  rho(xi_rho, xi_z)
#   <name>_density.png  ->  2*pi*xi_rho * rho(xi_rho, xi_z)
#
# Color map: black = 0, red = half of max, yellow = max.
# ---------------------------------------------------------------------------

if (!exists("output_dir")) output_dir = "output_cyl"
if (!exists("plot_dir")) plot_dir = "cyl_plots"
if (!exists("grid_type")) grid_type = "cyl"
system(sprintf("mkdir -p %s", plot_dir))
datafile = output_dir . "/dens.dat"
gridfile = "dens_grid.tmp"

# --- Build a pm3d-friendly copy with one blank-line-separated scan per block
if (grid_type eq "cyl") {
    system(sprintf("awk '/^[[:space:]]*#/{next} NF==0{next} \
    { for(i=1;i<=NF;i++) if ($i !~ /[EeDd]/ && $i ~ /[0-9][+-][0-9][0-9][0-9]$/) \
        $i=substr($i,1,length($i)-4) \"E\" substr($i,length($i)-3); \
      if (prev!=\"\" && $1!=prev) print \"\"; print; prev=$1 }' %s > %s", \
    datafile, gridfile))
} else {
    system(sprintf("awk '/^[[:space:]]*#/{next} NF==0{next} \
    { for(i=1;i<=NF;i++) if ($i !~ /[EeDd]/ && $i ~ /[0-9][+-][0-9][0-9][0-9]$/) \
        $i=substr($i,1,length($i)-4) \"E\" substr($i,length($i)-3); \
      r=sqrt($1*$1+$2*$2); if (prev!=\"\" && r<prev) print \"\"; print; prev=r }' %s > %s", \
    datafile, gridfile))
}

# --- Discover number of columns from the first data line
ncol = 0 + system(sprintf("awk '!/^[[:space:]]*#/ && NF>0 {print NF; exit}' %s", datafile))

# --- Terminal / style
cw  = 1800.0             # canvas width (px)
ch0 = 1500.0            # canvas height (px) used only for the geometry calibration below
set terminal pngcairo enhanced size cw,ch0 font ",28"
set pm3d map
set palette defined (0 "black", 0.5 "red", 1 "yellow")
set object 1 rectangle from graph 0,0 to graph 1,1 behind fillcolor rgb "black" fillstyle solid 1.0 noborder
set contour base                 # draw contour lines on the base plane
set cntrparam levels 10          # number of contour levels
unset clabel                     # no numeric labels on the contour lines
set tics out             # ticks point outward
set tics scale 0.5       # half the default tick length
set xtics offset 0,0.5   # move the x tic labels ~twice closer to the plot box
set xlabel "{/Symbol x}_{/Symbol r}" offset 0,0.5   # shift the x-axis name up by the same amount
set ylabel "{/Symbol x}_z"
set cblabel ""
set size ratio -1            # equal aspect for the two spatial axes
set autoscale fix

# --- Tighten the vertical white space without resizing the plot box.
#     "set size ratio -1" picks the natural (aspect-correct) box and pads it
#     with equal top/bottom white margins.  We render one throwaway frame to
#     read that box geometry back, then replace ratio -1 with fixed screen
#     margins that keep the box the very same size while halving the upper
#     white margin and trimming the lower one by "shift" (the x-label up-shift).
#     The canvas height is shrunk to match so the removed space is not re-added.
shift = 25.0             # px, matches the 0.5-character upward shift of the x labels
calfile = gridfile . ".cal.png"
stats gridfile using 3 nooutput              # representative color range ...
set cbrange [0:STATS_max]
set title "{/Symbol r}({/Symbol x}_{/Symbol r}, {/Symbol x}_z)"   # ... and a one-line title, so the
set output calfile                           # calibrated box matches the real frames
splot gridfile using 1:2:3 with pm3d notitle
unset output
boxL = GPVAL_TERM_XMIN
boxR = GPVAL_TERM_XMAX
topW = GPVAL_TERM_YMAX             # current upper white margin (px)
botW = ch0 - GPVAL_TERM_YMIN       # current lower white margin (px)
boxH = GPVAL_TERM_YMIN - GPVAL_TERM_YMAX
ch   = boxH + topW/2.0 + (botW - shift)   # new, shorter canvas height
unset size                  # box position is now fixed explicitly below
set terminal pngcairo enhanced size cw,ch font ",28"
set lmargin at screen boxL/cw
set rmargin at screen boxR/cw
set tmargin at screen (ch - topW/2.0)/ch
set bmargin at screen (botW - shift)/ch
system(sprintf("rm -f %s", calfile))

# --- Loop over the data columns (column 3 .. ncol)
do for [c=3:ncol] {
    name = system(sprintf("awk 'NR==1{gsub(/#/,\"\",$%d); print $%d}' %s", c, c, datafile))
    glabel = "{/Symbol r}_{" . substr(name, 4, strlen(name)) . "}"

    # The finite-mass nuclear density is unresolved beyond the origin on this
    # sample grid; the electronic density has appreciable support to about 4.
    if (name eq "rho1") {
        direct_rho = 0.0009
        direct_z = 0.0009
        axial_rho = 0.0010
        axial_z = 0.0010
        set format x "%.1te%T"
        set format y "%.1te%T"
        set xtics 0.0003
        set ytics 0.0003
        set pm3d interpolate 10,10
        unset contour
    } else {
        direct_rho = 1.0
        direct_z = 1.0
        axial_rho = 1.75
        axial_z = 1.75
        set format x "%g"
        set format y "%g"
        set xtics autofreq
        set ytics autofreq
        set pm3d interpolate 3,3
        set contour base
    }

    # ---- plot 1: rho itself ----------------------------------------------
    set autoscale x
    set autoscale y
    stats gridfile using c nooutput
    set cbrange [0:STATS_max]
    set xrange [0:direct_rho]
    set yrange [-direct_z:direct_z]
    set title sprintf("%s({/Symbol x}_{/Symbol r}, {/Symbol x}_z)", glabel)
    set output sprintf("%s/%s.png", plot_dir, name)
    splot gridfile using 1:2:c with pm3d notitle, \
          gridfile using 1:2:c with lines lc rgb "black" nosurface notitle
    unset output

    # ---- plot 2: probability density 2*pi*xi_rho*rho ---------------------
    set autoscale x
    set autoscale y
    stats gridfile using (2*pi*column(1)*column(c)) nooutput
    if (STATS_max > 0) {
        set cbrange [0:STATS_max]
        set xrange [0:axial_rho]
        set yrange [-axial_z:axial_z]
        set title sprintf("2{/Symbol p}{/Symbol x}_{/Symbol r} %s({/Symbol x}_{/Symbol r}, {/Symbol x}_z)", glabel)
        set output sprintf("%s/%s_axial_dist.png", plot_dir, name)
        splot gridfile using 1:2:(2*pi*column(1)*column(c)) with pm3d notitle, \
              gridfile using 1:2:(2*pi*column(1)*column(c)) with lines lc rgb "black" nosurface notitle
        unset output
    } else {
        print sprintf("Skipping %s_axial_dist.png (function is identically zero)", name)
    }
}

system(sprintf("rm -f %s", gridfile))
