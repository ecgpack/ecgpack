#!/usr/bin/env python3
"""Generate piecewise radial grids for density and correlation functions."""

from pathlib import Path

OUTPUT_DIRECTORY = Path(__file__).resolve().parent


class Grid:
    def __init__(self, counts, powers, x_values, gridtype):
        # Each grid must own its parameters so that refinement of one grid does
        # not modify the other grids through shared mutable objects.
        self.counts = list(counts)
        self.powers = list(powers)
        self.x_values = list(x_values)
        self.type = gridtype

    def points(self):
        """Yield points from each segment without duplicating shared boundaries."""
        for segment, (count, power) in enumerate(zip(self.counts, self.powers)):
            start = 0 if segment == 0 else 1
            x_min = self.x_values[segment]
            x_max = self.x_values[segment + 1]

            for index in range(start, count + 1):
                fraction = index / count
                yield x_min + (x_max - x_min) * fraction**power

    def add_tail(self, tail_count, tail_power, tail_end):
        """Add a tail segment to the grid."""
        self.counts.append(tail_count)
        self.powers.append(tail_power)
        self.x_values.append(tail_end)

    def increase_counts(self, num):
        """Multiply the number of intervals in each segment by num."""
        self.counts = [count * num for count in self.counts]

    def change_type(self, gridtype):
        """Change the grid type."""
        self.type = gridtype


def write_grid(path, grid):
    i=0
    with path.open("w", encoding="utf-8") as grid_file:
        for point in grid.points():
            grid_file.write(f"{point:.12f}\n")
            i+=1
    return i

def main():
    counts = [200, 250, 200, 200, 200, 200]
    powers = [1.0, 1.0, 1.1, 1.0, 1.0, 1.0]
    x_values = [0.0, 0.002, 1.0, 10.0, 20.0, 48.0,
        64.0]

    sample_grid = Grid(counts, powers, x_values, "sample_radial")
    path = OUTPUT_DIRECTORY / f"{sample_grid.type}_ecg_grid.dat"
    point_count = write_grid(path, sample_grid)
    print(f"Wrote {path.name:<19} ({point_count} points)")


if __name__ == "__main__":
    main()
