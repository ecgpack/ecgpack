#!/usr/bin/env python3
"""Generate high-precision mpmath reference values for K_0(x) and K_1(x)."""

from __future__ import annotations

import argparse
from pathlib import Path

import mpmath as mp


DEFAULT_DPS = 120
DEFAULT_OUTPUT_DIGITS = 90
DEFAULT_X_VALUES = (
    "1e-40",
    "1e-32",
    "1e-16",
    "1e-8",
    "1e-6",
    "1e-3",
    "1e-1",
    "5e-1",
    "1",
    "2",
    "5",
    "10",
    "20",
    "40",
    "100",
    "1e+4",
    "1e+6",
    "1e+10",
    "1e+16",
    "1e+32",
)
DEFAULT_OUTPUT = Path(__file__).with_name("bessel_k01_mpmath_reference.dat")


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate a table containing K_0(x), K_1(x), exp(x)*K_0(x), "
            "and exp(x)*K_1(x) using mpmath."
        )
    )
    parser.add_argument(
        "--dps",
        type=int,
        default=DEFAULT_DPS,
        help=f"mpmath working precision in decimal digits (default: {DEFAULT_DPS})",
    )
    parser.add_argument(
        "--digits",
        type=int,
        default=DEFAULT_OUTPUT_DIGITS,
        help=(
            "significant decimal digits written for each value "
            f"(default: {DEFAULT_OUTPUT_DIGITS})"
        ),
    )
    parser.add_argument(
        "--x-values",
        nargs="+",
        default=DEFAULT_X_VALUES,
        metavar="X",
        help="positive x values; decimal strings are read without binary-float conversion",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"output table path (default: {DEFAULT_OUTPUT})",
    )
    return parser.parse_args()


def validate_arguments(arguments: argparse.Namespace) -> None:
    if arguments.dps < 20:
        raise ValueError("--dps must be at least 20")
    if arguments.digits < 1:
        raise ValueError("--digits must be positive")
    if arguments.digits > arguments.dps - 10:
        raise ValueError(
            "--dps must exceed --digits by at least 10 guard digits"
        )
    if not arguments.x_values:
        raise ValueError("at least one x value is required")


def format_value(value: mp.mpf, digits: int) -> str:
    return mp.nstr(value, digits, min_fixed=0, max_fixed=0)


def generate_table(
    x_strings: list[str] | tuple[str, ...], dps: int, output_digits: int
) -> str:
    rows = [
        f"# Modified Bessel reference values generated with mpmath {mp.__version__}.",
        f"# mp.dps = {dps}; values below retain {output_digits} significant decimal digits.",
        "# Columns: x  K_0(x)  K_1(x)  exp(x)*K_0(x)  exp(x)*K_1(x)",
    ]

    with mp.workdps(dps):
        for x_string in x_strings:
            x = mp.mpf(x_string)
            if x <= 0:
                raise ValueError(f"x must be positive, received {x_string!r}")

            k0 = mp.besselk(0, x)
            k1 = mp.besselk(1, x)
            scale = mp.exp(x)
            values = (x, k0, k1, scale * k0, scale * k1)
            rows.append(
                "  ".join(format_value(value, output_digits) for value in values)
            )

    return "\n".join(rows) + "\n"


def main() -> None:
    arguments = parse_arguments()
    validate_arguments(arguments)
    table = generate_table(arguments.x_values, arguments.dps, arguments.digits)

    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(table, encoding="utf-8")
    print(f"Wrote {len(arguments.x_values)} reference rows to {arguments.output}")


if __name__ == "__main__":
    main()
