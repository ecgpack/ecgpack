#!/usr/bin/env python3
"""Generate a piecewise cylindrical grid for RG_2D densities."""

from __future__ import annotations

from pathlib import Path


OUTPUT_DIRECTORY = Path(__file__).resolve().parent


class CylindricalGrid:
    def __init__(self, rho_counts, rho_powers, rho_values,
                 z_counts, z_powers, z_values, gridtype):
        self.rho_counts = list(rho_counts)
        self.rho_powers = list(rho_powers)
        self.rho_values = list(rho_values)
        self.z_counts = list(z_counts)
        self.z_powers = list(z_powers)
        self.z_values = list(z_values)
        self.type = gridtype

    @staticmethod
    def increasing_axis_points(counts, powers, values):
        """Return one nonnegative axis without duplicate boundaries."""
        points = []
        for segment, (count, power) in enumerate(zip(counts, powers)):
            start = 0 if segment == 0 else 1
            x_min = values[segment]
            x_max = values[segment + 1]
            for index in range(start, count + 1):
                fraction = index / count
                points.append(x_min + (x_max - x_min) * fraction**power)
        return points

    def rho_points(self):
        """Return the nonnegative cylindrical-radius axis."""
        return self.increasing_axis_points(
            self.rho_counts, self.rho_powers, self.rho_values
        )

    def z_points(self):
        """Return an increasing axis from negative to positive z."""
        positive = self.increasing_axis_points(
            self.z_counts, self.z_powers, self.z_values
        )
        return [-point for point in reversed(positive[1:])] + positive

    def points(self):
        """Yield the tensor product in fixed-rho scans."""
        rhos = self.rho_points()
        zs = self.z_points()
        for rho in rhos:
            for z in zs:
                yield rho, z


def write_grid(path, grid):
    point_count = 0
    with path.open("w", encoding="utf-8") as grid_file:
        for rho, z in grid.points():
            grid_file.write(f"{rho:.12f} {z:.12f}\n")
            point_count += 1
    return point_count


def main():
    axis_values = [0.0, 0.002, 1.0, 8.0, 24.0, 48.0]
    axis_powers = [1.0, 1.1, 1.0, 1.0, 1.0]
    rho_counts = [42, 22, 15, 8, 4]
    z_counts = [42, 22, 15, 8, 4]

    sample_grid = CylindricalGrid(
        rho_counts, axis_powers, axis_values,
        z_counts, axis_powers, axis_values,
        "sample_cyl",
    )
    path = OUTPUT_DIRECTORY / f"{sample_grid.type}_ecg_grid.dat"
    point_count = write_grid(path, sample_grid)
    print(f"Wrote {path.name} ({point_count} points)")


if __name__ == "__main__":
    main()
