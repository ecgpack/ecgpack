#!/usr/bin/env python3
"""Generate a piecewise radial-and-angle grid for RG_1P densities."""

from __future__ import annotations

import math
from pathlib import Path


OUTPUT_DIRECTORY = Path(__file__).resolve().parent


class RadialAngularGrid:
    def __init__(self, radial_counts, radial_powers, radial_values,
                 phi_count, gridtype):
        self.radial_counts = list(radial_counts)
        self.radial_powers = list(radial_powers)
        self.radial_values = list(radial_values)
        self.phi_count = phi_count
        self.type = gridtype

    @staticmethod
    def increasing_axis_points(counts, powers, values):
        """Yield one nonnegative axis without duplicate segment boundaries."""
        for segment, (count, power) in enumerate(zip(counts, powers)):
            start = 0 if segment == 0 else 1
            x_min = values[segment]
            x_max = values[segment + 1]

            for index in range(start, count + 1):
                fraction = index / count
                yield x_min + (x_max - x_min) * fraction**power

    def radial_points(self):
        """Yield radii."""
        yield from self.increasing_axis_points(
            self.radial_counts, self.radial_powers, self.radial_values
        )

    def angular_points(self):
        """Yield angles in degrees."""
        yield from self.increasing_axis_points(
            [self.phi_count], [1], [0, 180]
        )

    def points(self):
        """Yield radii, angles, and the corresponding cylindrical points."""
        radii = list(self.radial_points())
        phis = list(self.angular_points())
        for phi in phis:
            angle = math.radians(phi)
            for radius in radii:
                yield (
                    radius,
                    phi,
                    radius * math.sin(angle),
                    radius * math.cos(angle),
                )

    def increase_rad_counts(self, factor: int):
        self.radial_counts = [count * factor for count in self.radial_counts]

    def increase_phi_count(self, factor: int): self.phi_count *= factor


def write_grids(ecg_path, radial_path, grid):
    point_count = 0
    with (
        ecg_path.open("w", encoding="utf-8") as ecg_file,
        radial_path.open("w", encoding="utf-8") as radial_file,
    ):
        radial_file.write("#r phi_deg grid_rho grid_z\n")
        for radius, phi, p_rho, p_z in grid.points():
            ecg_file.write(f"{p_rho:.12f} {p_z:.12f}\n")
            radial_file.write(
                f"{radius:.12f} {phi:.8f} {p_rho:.12f} {p_z:.12f}\n"
            )
            point_count += 1
    return point_count


def write_files(grid: RadialAngularGrid):
    """Write the ECGPACK grid and an angle-labelled reference table."""
    ecg_path = OUTPUT_DIRECTORY / f"{grid.type}_ecg_grid.dat"
    radial_path = OUTPUT_DIRECTORY / f"{grid.type}_coords.dat"

    count = write_grids(ecg_path, radial_path, grid)
    print(f"Wrote {ecg_path.name} ({count} points)")
    print(f"Wrote {radial_path.name} ({count} points)")


def main():
    radial_values = [
        0.0, 0.002, 1.0, 8.0, 24.0, 48.0
    ]
    radial_counts = [
        100, 100, 100, 50, 50
    ]
    radial_powers = [1, 1.1, 1, 1.1, 1]
    phi_count = 60

    sample_grid = RadialAngularGrid(
        radial_counts, radial_powers, radial_values, phi_count, "sample_radial"
    )
    write_files(sample_grid)


if __name__ == "__main__":
    main()
