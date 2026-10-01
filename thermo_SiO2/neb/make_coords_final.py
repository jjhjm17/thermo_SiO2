"""Create a LAMMPS ``coords.final`` coordinate file from a LAMMPS data file.

LAMMPS NEB accepts a comment line, the atom count, then ``id x y z`` records.
The source data file may also contain atom types and image flags; those fields
are deliberately omitted here while preserving the atom IDs and coordinates.
chat gpt
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_INPUT = SCRIPT_DIR / '../c.final/relaxed.dataf'
DEFAULT_OUTPUT = SCRIPT_DIR / 'coords.final'
ATOM_COUNT_RE = re.compile(r'^\s*(\d+)\s+atoms\s*$', re.IGNORECASE)
ATOMS_SECTION_RE = re.compile(r'^\s*Atoms(?:\s|#|$)', re.IGNORECASE)


def read_atom_count(lines: list[str], source: Path) -> int:
    """Return the declared atom count from a LAMMPS data-file header."""
    for line in lines:
        match = ATOM_COUNT_RE.match(line)
        if match:
            return int(match.group(1))
    raise ValueError(f'{source}: no "<count> atoms" header found')


def read_atom_coordinates(source: Path) -> list[tuple[int, str, str, str]]:
    """Extract exactly the declared records from the ``Atoms`` section."""
    lines = source.read_text(encoding='utf-8').splitlines()
    atom_count = read_atom_count(lines, source)

    try:
        atoms_start = next(
            index + 1 for index, line in enumerate(lines)
            if ATOMS_SECTION_RE.match(line)
        )
    except StopIteration as error:
        raise ValueError(f'{source}: no Atoms section found') from error

    coordinates: list[tuple[int, str, str, str]] = []
    atom_ids: set[int] = set()
    for line_number, line in enumerate(lines[atoms_start:], start=atoms_start + 1):
        if not line.strip():
            continue
        fields = line.split()
        if len(fields) < 5:
            raise ValueError(
                f'{source}:{line_number}: expected "id type x y z" in Atoms section',
            )
        try:
            atom_id = int(fields[0])
            float(fields[2])
            float(fields[3])
            float(fields[4])
        except ValueError as error:
            raise ValueError(
                f'{source}:{line_number}: invalid atom ID or coordinate',
            ) from error
        if atom_id in atom_ids:
            raise ValueError(f'{source}:{line_number}: duplicate atom ID {atom_id}')

        atom_ids.add(atom_id)
        coordinates.append((atom_id, fields[2], fields[3], fields[4]))
        if len(coordinates) == atom_count:
            return coordinates

    raise ValueError(
        f'{source}: Atoms section contains {len(coordinates)} records; '
        f'header declares {atom_count}',
    )


def write_neb(coordinates: list[tuple[int, str, str, str]], destination: Path) -> None:
    """Write coordinates in LAMMPS NEB final-coordinate format."""
    with destination.open('w', encoding='utf-8') as file:
        file.write('# Generated from LAMMPS data file\n\n')
        file.write(f'{len(coordinates)}\n')
        for atom_id, x, y, z in coordinates:
            file.write(f'{atom_id} {x} {y} {z}\n')


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', nargs='?', type=Path, default=DEFAULT_INPUT,
                        help=f'Lammps data file (default: {DEFAULT_INPUT})')
    # positional arguments (input) cf. optional: '--input'
    parser.add_argument('output', nargs='?', type=Path, default=DEFAULT_OUTPUT,
                        help=f'NEB coordinate file (default: {DEFAULT_OUTPUT})')
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    coordinates = read_atom_coordinates(args.input)
    write_neb(coordinates, args.output)
    print(f'Wrote {len(coordinates)} coordinates to {args.output}')


if __name__ == '__main__':
    main()
