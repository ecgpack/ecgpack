#!/usr/bin/env gnuplot
# ---------------------------------------------------------------------------
# plot_cf.gp
#
# Reads cf.dat:
#   col 1 = xi_rho   (x-axis)
#   col 2 = xi_z     (y-axis)
#   col 3..N = pair correlation functions g_i (names taken from the header)
#
# For every data column it produces two PNG contour/density plots:
#   <name>.png          ->  g(xi_rho, xi_z)
#   <name>_density.png  ->  2*pi*xi_rho * g(xi_rho, xi_z)
#
# Color map: black = 0, red = half of max, yellow = max.
# ---------------------------------------------------------------------------

datafile = "cf.dat"
gridfile = "cf_grid.tmp"

# --- Build a pm3d-friendly copy: drop comments, blank line between xi_rho blocks
system(sprintf("awk '/^[[:space:]]*#/{next} NF==0{next} \
{ if (prev!=\"\" && $1!=prev) print \"\"; print; prev=$1 }' %s > %s", \
datafile, gridfile))

# --- Discover number of columns from the first data line
ncol = 0 + system(sprintf("awk '!/^[[:space:]]*#/ && NF>0 {print NF; exit}' %s", datafile))

# --- Terminal / style
cw  = 1800.0             # canvas width (px)
ch0 = 1500.0            # canvas height (px) used only for the geometry calibration below
set terminal pngcairo enhanced size cw,ch0 font ",28"
set pm3d map
set palette defined (0 "black", 0.5 "red", 1 "yellow")
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
set title "g({/Symbol x}_{/Symbol r}, {/Symbol x}_z)"   # ... and a one-line title, so the
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
    glabel = "g_{" . substr(name, 2, strlen(name)) . "}"

    # ---- plot 1: g itself ------------------------------------------------
    stats gridfile using c nooutput
    set cbrange [0:STATS_max]
    set title sprintf("%s({/Symbol x}_{/Symbol r}, {/Symbol x}_z)", glabel)
    set output sprintf("%s.png", name)
    splot gridfile using 1:2:c with pm3d notitle, \
          gridfile using 1:2:c with lines lc rgb "black" nosurface notitle
    unset output

    # ---- plot 2: probability density 2*pi*xi_rho*g -----------------------
    stats gridfile using (2*pi*column(1)*column(c)) nooutput
    if (STATS_max > 0) {
        set cbrange [0:STATS_max]
        set title sprintf("2{/Symbol p}{/Symbol x}_{/Symbol r} %s({/Symbol x}_{/Symbol r}, {/Symbol x}_z)", glabel)
        set output sprintf("%s_axial_dist.png", name)
        splot gridfile using 1:2:(2*pi*column(1)*column(c)) with pm3d notitle, \
              gridfile using 1:2:(2*pi*column(1)*column(c)) with lines lc rgb "black" nosurface notitle
        unset output
    } else {
        print sprintf("Skipping %s_axial_dist.png (function is identically zero)", name)
    }
}

system(sprintf("rm -f %s", gridfile))
