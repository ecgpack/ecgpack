#!/usr/bin/env bash
#
# Generates two radial grid files (one r value per row):
#   cf_grid.dat    - uniform 2D grid from xi_rho_min_cf to xi_rho_max_cf with step xi_rho_step_cf (note that $xi_rho=\sqrt{(xi_x)^2+(xi_y)^2}$ here)
#                    and z from xi_z_min_cf to xi_z_max_cf with step xi_z_step_cf
#   dens_grid.dat  - uniform 2D grid from xi_rho_min_dens to xi_rho_max_dens with step xi_rho_step_dens (note that $xi_rho=\sqrt{(xi_x)^2+(xi_y)^2}$ here)
#                    and z from xi_z_min_dens to xi_z_max_dens with step xi_z_step_dens

set -euo pipefail

# ---- cf grid parameters ----
xi_rho_min_cf=0.0       # minimum xi_rho (must be non-negative)
xi_rho_max_cf=2.0       # maximum xi_rho
xi_rho_step_cf=0.02     # step size of the grid
xi_z_min_cf=-1.5        # minimum xi_z
xi_z_max_cf=1.5         # maximum xi_z
xi_z_step_cf=0.02       # step size of the grid

# ---- dens grid parameters ----
xi_rho_min_dens=0.0       # minimum xi_rho (must be non-negative)
xi_rho_max_dens=2.0       # maximum xi_rho
xi_rho_step_dens=0.02     # step size of the grid
xi_z_min_dens=-1.5        # minimum xi_z
xi_z_max_dens=1.5       # maximum xi_z
xi_z_step_dens=0.02       # step size of the grid

# ---- output format ----
# The leading + makes positive and negative numbers have the same width,
# e.g. +1.500000 and -1.500000.
fmt="%+.6f"

# ---- cf_grid.dat: ----
for r in $(seq -f "$fmt" "$xi_rho_min_cf" "$xi_rho_step_cf" "$xi_rho_max_cf"); do
  for z in $(seq -f "$fmt" "$xi_z_min_cf" "$xi_z_step_cf" "$xi_z_max_cf"); do
    echo "$r $z"
  done
done > cf_grid.dat

# ---- dens_grid.dat: ----
for r in $(seq -f "$fmt" "$xi_rho_min_dens" "$xi_rho_step_dens" "$xi_rho_max_dens"); do
  for z in $(seq -f "$fmt" "$xi_z_min_dens" "$xi_z_step_dens" "$xi_z_max_dens"); do
    echo "$r $z"
  done
done > dens_grid.dat

echo "Wrote cf_grid.dat   ($(wc -l < cf_grid.dat) points)"
echo "Wrote dens_grid.dat ($(wc -l < dens_grid.dat) points)"