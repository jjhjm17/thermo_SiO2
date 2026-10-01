"""Collect the final geometry of each LAMMPS NEB image as a VASP POSCAR.

By default, input files are read from (for example) ``../d.neb/relax.dump.*`` and output
files are written here as ``POSCAR.<image-number>``.  The dump suffixes are
sorted numerically, so (for example) image 10 follows image 9 rather than
image 1.

change  cell to (for example) [[19, 0, 0], [0, 16, 0], [0, 0, 16]] so that the vacuum size is at least 10 Ang.


"""

from __future__ import annotations

import argparse
from pathlib import Path

from ase.io import read, write
# from thermo_SiO2.mlip.write_config import write_cfg_SiO2


SPECORDER = ("Si", "O", "H", "Al")
# CELL = [[19, 0, 0], [0, 16, 0], [0, 0, 16]] 
NEB_DIR = '../d.e.neb/'


def image_number(path: Path) -> int:
    """Return the numeric suffix of a ``relax.dump.<image>`` file."""
    try:
        return int(path.name.removeprefix("relax.dump."))
    except ValueError as error:
        raise ValueError(f"Expected a numeric image suffix: {path}") from error


def collect_images(input_dir: Path, output_dir: Path) -> list[Path]:
    """Write the last snapshot of every NEB dump and return the outputs."""
    dumps = sorted(input_dir.glob("relax.dump.*"), key=image_number)
    if not dumps:
        raise FileNotFoundError(f"No files matching {input_dir / 'relax.dump.*'}")

    output_dir.mkdir(parents=True, exist_ok=True)
    outputs = []
    cfgs = []
    for dump in dumps:
        image = image_number(dump)
        atoms = read(dump, index=-1, format="lammps-dump-text", specorder=SPECORDER)

        # LAMMPS type IDs map to SPECORDER; arrange atoms in contiguous species
        # blocks, as expected by a conventional VASP POSCAR.
        order = sorted(
            range(len(atoms)), key=lambda index: SPECORDER.index(atoms[index].symbol)
        )
        atoms = atoms[order]
        # atoms.set_cell(CELL)

        # output = output_dir / f"POSCAR.{image}"
        # write(output, atoms, format="vasp", direct=True, vasp5=True, sort=False)
        cfgs.append(atoms)
        # outputs.append(output)

    # write_cfg_SiO2('neb.cfg', cfgs, atom_symbols='Si O H Al')
    write('neb.xyz', cfgs)
    outputs = ['neb.xyz']
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "input_dir",
        nargs="?",
        type=Path,
        default=Path(NEB_DIR),
        help="directory containing relax.dump.* files (default: ../d.neb)",
    )
    parser.add_argument(
        "output_dir",
        nargs="?",
        type=Path,
        default=Path("."),
        help="directory for POSCAR.<image> files (default: current directory)",
    )
    args = parser.parse_args()

    for output in collect_images(args.input_dir, args.output_dir):
        print(output)


if __name__ == "__main__":
    main()
