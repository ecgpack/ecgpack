#!/usr/bin/env bash
#
# Generates two radial grid files (one r value per row):
#   cf_grid.dat   - uniform grid from r_min to r_max with step r_step
#   dens_grid.dat - fine grid from r_min to r_nuc_max (step r_nuc_step),
#                  then a coarser grid up to r_max (step r_step).
#
# The fine inner region is meant to resolve the density of the first
# particle, which is expected to be highly localized around r = 0.

set -euo pipefail

# ---- coarse grid parameters ----
r_min=0.0          # minimum r
r_max=9.0          # maximum r
r_step=0.02        # step size of the coarse (outer) grid

# ---- fine grid parameters ----
r_nuc_max=0.0009      # upper bound of the fine (inner) region
r_nuc_step=0.000005   # step size of the fine (inner) grid

# ---- output format ----
fmt="%.6f"

# ---- cfgrid.dat: single uniform grid ----
seq -f "$fmt" "$r_min" "$r_step" "$r_max" > cf_grid.dat

# ---- densgrid.dat: fine inner grid + coarse outer grid ----
# Inner: r_min .. r_nuc_max with step r_nuc_step
seq -f "$fmt" "$r_min" "$r_nuc_step" "$r_nuc_max" > dens_grid.dat
# Outer: continue from (r_nuc_max + r_step) .. r_max with step r_step
#        so the r_nuc_max point is not duplicated
r_outer_start=$(awk -v a="$r_nuc_max" -v s="$r_step" 'BEGIN{printf "%.10g", a + s}')
seq -f "$fmt" "$r_outer_start" "$r_step" "$r_max" >> dens_grid.dat

echo "Wrote cf_grid.dat   ($(wc -l < cf_grid.dat) points)"
echo "Wrote dens_grid.dat ($(wc -l < dens_grid.dat) points)"
