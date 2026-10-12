#!/usr/bin/env python3
"""Generate the radial grid for the nucleus-nucleus correlation function."""

from pathlib import Path

OUTPUT_DIRECTORY = Path(__file__).resolve().parent


class Grid:
    def __init__(self, counts, powers, x_values, name):
        # Each segment i runs from x_values[i] to x_values[i+1] with counts[i]
        # intervals; power 1.0 gives equal steps, power > 1.0 crowds the points
        # towards the start of the segment.
        self.counts = list(counts)
        self.powers = list(powers)
        self.x_values = list(x_values)
        self.name = name

    def points(self):
        """Yield points from each segment without duplicating shared boundaries."""
        for segment, (count, power) in enumerate(zip(self.counts, self.powers)):
            start = 0 if segment == 0 else 1
            x_min = self.x_values[segment]
            x_max = self.x_values[segment + 1]

            for index in range(start, count + 1):
                fraction = index / count
                yield x_min + (x_max - x_min) * fraction**power


def write_grid(path, grid):
    i = 0
    with path.open("w", encoding="utf-8") as grid_file:
        for point in grid.points():
            grid_file.write(f"{point:.12f}\n")
            i += 1
    return i


def main():
    # One segment from 0 to 6 bohr with 600 equal steps (601 points, step 0.01).
    # The internuclear distance of HD is about 1.4 bohr; the grid reaches far
    # enough into the tail for the integrals of g1 to be converged.
    counts = [600]
    powers = [1.0]
    x_values = [0.0, 6.0]

    grid = Grid(counts, powers, x_values, "cf_grid")
    path = OUTPUT_DIRECTORY / f"{grid.name}.dat"
    point_count = write_grid(path, grid)
    print(f"Wrote {path.name} ({point_count} points)")


if __name__ == "__main__":
    main()
