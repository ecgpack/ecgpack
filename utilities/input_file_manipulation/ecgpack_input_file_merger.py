#!/bin/python3
# This Python script merges two ECGPACK input files into one. It takes a
# (sub)range of basis functions from each of the two input files and writes a
# new input file whose basis is the concatenation of the two selected ranges.
#
# Usage:
#
#   ecgpack_input_file_merger.py -file1 <file1name> -range1 <X1>-<Y1> \
#                                       -file2 <file2name> -range2 <X2>-<Y2> \
#                                       -file3 <file3name>
#
#   <file1name>, <file2name> - names of the two input files that are merged
#   <file3name>              - name of the resulting (merged) input file
#   <X1>-<Y1>, <X2>-<Y2>     - ranges of basis functions taken from each file
#                              (two integers separated by a dash, e.g. 5-20)
#
# Arguments -file1, -file2, and -file3 are required. Arguments -range1 and
# -range2 are optional; if a range is not supplied it defaults to the full
# basis of the corresponding file (functions 1 through K, where K is the basis
# size of that file).
#
# Examples:
#
#   ecgpack_input_file_merger.py -file1 a.txt -file2 b.txt -file3 c.txt
#   ecgpack_input_file_merger.py -file1 a.txt -range1 5-20 -file2 b.txt -range2 1-8 -file3 c.txt
#
# In the resulting file:
#   * the header keywords and their values are copied from <file1name>, except
#     for BASIS_SIZE, which is updated to the size of the merged basis;
#   * the command list contains a single command, "EXPC_VALS G <mergedbasissize>";
#   * the history is emptied (all energies set to zero and all counters set to
#     zero), except for its last line, which carries the merged basis size and
#     the current energy taken from the last history line of <file1name>;
#   * the basis functions are the selected range from <file1name> followed by
#     the selected range from <file2name>, renumbered from 1 to the merged size.
#
# The script verifies that the two input files describe the same number of
# particles (keyword PARTICLES) and that their basis-function lines have the
# same number of entries (which guards against accidentally merging files that
# use different basis types). Real numbers are copied verbatim as strings, so
# their length/precision (double, extended, quadruple) is preserved.

import argparse
import re
import sys


def error(message):
    """Print an error message and terminate the script."""
    sys.stderr.write("Error: " + message + "\n")
    sys.exit(1)


def parse_file(path):
    """Split an ECGPACK input file into its five parts.

    Returns (header, excerpt, commands, history, basis, ruler), where each of
    the first five is a list of lines (without trailing newlines) and ruler is
    the exact ruler-line string used in the file.
    """
    try:
        with open(path) as f:
            lines = [ln.rstrip("\n") for ln in f]
    except OSError as exc:
        error("cannot open file '{}': {}".format(path, exc))

    # Ruler lines are those consisting solely of '=' characters (surrounding
    # whitespace ignored). A complete input file contains exactly four of them.
    rulers = [i for i, ln in enumerate(lines)
              if ln.strip() and set(ln.strip()) == {"="}]
    if len(rulers) != 4:
        error("file '{}' does not look like a complete input file "
              "(expected 4 ruler lines, found {}).".format(path, len(rulers)))

    header = lines[:rulers[0]]
    excerpt = lines[rulers[0] + 1:rulers[1]]
    commands = lines[rulers[1] + 1:rulers[2]]
    history = lines[rulers[2] + 1:rulers[3]]
    basis = lines[rulers[3] + 1:]

    # Discard any trailing blank lines in the basis-function section.
    while basis and basis[-1].strip() == "":
        basis.pop()

    return header, excerpt, commands, history, basis, lines[rulers[0]]


def get_header_int(header, keyword, path):
    """Return the integer value that follows a header keyword."""
    for ln in header:
        toks = ln.split()
        if toks and toks[0] == keyword:
            try:
                return int(toks[1])
            except (IndexError, ValueError):
                error("cannot read the value of '{}' in file '{}'."
                      .format(keyword, path))
    error("keyword '{}' not found in the header of file '{}'."
          .format(keyword, path))


def basis_entry_count(basis, path):
    """Return the (uniform) number of whitespace-separated entries per line."""
    counts = set(len(ln.split()) for ln in basis if ln.strip())
    if not counts:
        error("file '{}' contains no basis functions to merge.".format(path))
    if len(counts) != 1:
        error("basis-function lines in file '{}' do not all have the same "
              "number of entries.".format(path))
    return counts.pop()


def parse_range(range_str, ksize, which):
    """Parse a 'X-Y' range string and validate it against the basis size."""
    m = re.match(r"^\s*(\d+)\s*-\s*(\d+)\s*$", range_str)
    if not m:
        error("range '{}' for {} is not of the form X-Y (two integers "
              "separated by a dash).".format(range_str, which))
    x, y = int(m.group(1)), int(m.group(2))
    if x < 1 or y < 1:
        error("range '{}' for {} must use positive function numbers."
              .format(range_str, which))
    if x > y:
        error("range '{}' for {} has its start greater than its end."
              .format(range_str, which))
    if y > ksize:
        error("range '{}' for {} exceeds the basis size ({}) of that file."
              .format(range_str, which, ksize))
    return x, y


def history_line(size, energy_str):
    """Format a single history line: size, energy, and three zero counters.

    The energy string is written verbatim (its precision/length is preserved).
    A sign slot is kept so that positive numbers align with negative ones.
    """
    if energy_str.startswith("-") or energy_str.startswith("+"):
        efield = energy_str
    else:
        efield = " " + energy_str
    return "{:>7} {}{:>7}{:>7}{:>7}".format(size, efield, 0, 0, 0)


def renumber(line, new_index):
    """Replace the leading integer index of a basis-function line."""
    m = re.match(r"\s*\d+(.*)$", line)
    return "{:>7}".format(new_index) + m.group(1)


def main():
    parser = argparse.ArgumentParser(
        description="Merge two ECGPACK input files into one.")
    parser.add_argument("-file1", required=True, help="First input file.")
    parser.add_argument("-range1", default=None,
                        help="Range X1-Y1 of basis functions taken from file1 "
                             "(default: the whole basis).")
    parser.add_argument("-file2", required=True, help="Second input file.")
    parser.add_argument("-range2", default=None,
                        help="Range X2-Y2 of basis functions taken from file2 "
                             "(default: the whole basis).")
    parser.add_argument("-file3", required=True, help="Resulting merged file.")
    args = parser.parse_args()

    header1, _, _, history1, basis1, ruler1 = parse_file(args.file1)
    header2, _, _, history2, basis2, ruler2 = parse_file(args.file2)

    # The number of particles must be the same in both files.
    np1 = get_header_int(header1, "PARTICLES", args.file1)
    np2 = get_header_int(header2, "PARTICLES", args.file2)
    if np1 != np2:
        error("the number of particles differs between the two files "
              "({} in '{}' vs {} in '{}').".format(np1, args.file1,
                                                    np2, args.file2))

    # The basis-function lines must have the same number of entries; otherwise
    # the two files most likely use different basis types.
    nentry1 = basis_entry_count(basis1, args.file1)
    nentry2 = basis_entry_count(basis2, args.file2)
    if nentry1 != nentry2:
        error("basis-function lines have a different number of entries "
              "({} in '{}' vs {} in '{}'); the files probably use different "
              "basis types.".format(nentry1, args.file1, nentry2, args.file2))

    k1 = len(basis1)
    k2 = len(basis2)

    x1, y1 = parse_range(args.range1, k1, "file1") if args.range1 else (1, k1)
    x2, y2 = parse_range(args.range2, k2, "file2") if args.range2 else (1, k2)

    selected = basis1[x1 - 1:y1] + basis2[x2 - 1:y2]
    merged_size = len(selected)

    # Current energy for the merged file: the energy (second word) taken from
    # the last history line of file1.
    if not history1:
        error("file '{}' has an empty history section.".format(args.file1))
    try:
        energy1 = history1[-1].split()[1]
    except IndexError:
        error("cannot read the energy from the last history line of '{}'."
              .format(args.file1))

    # Assemble the merged file.
    out = []

    # Header: copy from file1, updating only BASIS_SIZE.
    for ln in header1:
        toks = ln.split()
        if toks and toks[0] == "BASIS_SIZE":
            m = re.match(r"(\s*BASIS_SIZE\s*)\d+.*$", ln)
            out.append(m.group(1) + str(merged_size))
        else:
            out.append(ln)

    # Header history excerpt (same form as the last history line).
    out.append(ruler1)
    out.append(history_line(merged_size, energy1))
    out.append(ruler1)

    # Command list: a single EXPC_VALS command.
    out.append(" EXPC_VALS G " + str(merged_size))

    # History: all zeros, except the last line which carries energy1.
    out.append(ruler1)
    for i in range(1, merged_size + 1):
        if i == merged_size:
            out.append(history_line(i, energy1))
        else:
            out.append(history_line(i, "0.0000000000000000E+00"))

    # Basis functions: file1 range first, then file2 range, renumbered 1..K.
    out.append(ruler1)
    for i, ln in enumerate(selected, start=1):
        out.append(renumber(ln, i))

    try:
        with open(args.file3, "w") as f:
            f.write("\n".join(out) + "\n")
    except OSError as exc:
        error("cannot write file '{}': {}".format(args.file3, exc))

    sys.stderr.write(
        "Merged {} function(s) [{}-{}] from '{}' and {} function(s) [{}-{}] "
        "from '{}' into '{}' ({} functions total).\n".format(
            y1 - x1 + 1, x1, y1, args.file1,
            y2 - x2 + 1, x2, y2, args.file2,
            args.file3, merged_size))


if __name__ == "__main__":
    main()
