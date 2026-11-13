"""
human doc
=========

This script reads the final energies of (for example) ../d.neb/log.lammps.* files and plot
Energy (eV) vs Reaction coordinate (Ang) to file fig_neb.pdf.

Also plot DFT energies in (for example) ../f.DFT/calc/*/OUTCAR folders, sorted.
(old) Read energy(sigma->0) (E0).
(new) Read TOTEN energy, free energy, not sigma->0, because mlip cfg has this energy.

How to find the last line?
After getting the header, check the last line and see if it is valid, that is, it has the same number of whitespace-
  separated fields.

--neb-lmp-dir   : for reading lammps output
--neb-grace-xyz: read xyz file for GRACE NEB, read free energy, plot. (Only one of --neb-lmp-dir or --neb-grace-xyz are needed.)
--dft-dir : DFT folder
--single-point-syn-mtp a.cfg : mlip cfg file of single point energies of syn-MTP.
            use thermo_SiO2.io.read_sil to get energies. 
--single-point-eq-mtp a.cfg : mlip cfg file of single point energies of eq-MTP.
--atom-symbols: used for --neb-mtp_single-point

plot_data.txt: round to 5 decimal digits (eV)

Save the ase Atoms objects to neb_DFT.xyz.

Keep comments in this comment block, these are written by a person.

Now the docstring by an agent starts:
=====================================

DFT energy: TOTEN of the last ionic step of each ``<dft-dir>/*/OUTCAR``, read by ASE
(format 'vasp-out') as ``free_energy`` (the force-consistent energy).
The OUTCAR folders are sorted by name (calc/0000, calc/0001, ...), and folder i is
compared with NEB image i. The number of OUTCARs must equal the number of NEB images.

--neb-grace-xyz: an extxyz file with all NEB images (e.g. neb_b_climb.xyz from run_neb.py).
  Energy per image: header ``free_energy=`` (read via ASE as the force-consistent energy).
  Reaction coordinate: cumulative image-to-image distance in 3N-dimensional space,
  using the minimum-image convention, normalized to [0, 1]. This is the same
  definition as ``print_neb_status`` in run_neb.py and LAMMPS's RDn.

--single-point-syn-mtp / --single-point-eq-mtp: MLIP cfg files with syn-MTP / eq-MTP
  single-point energies of the same NEB images (same order, same number of images).
  Read with thermo_SiO2.io.read_sil using --atom-symbols (default 'Si O H Al').
  Plotted on the same reaction coordinate as the NEB. Both are optional.

Every curve is plotted relative to its own first image.

Output: fig_neb.pdf and plot_data.txt. Columns: reaction coordinate, then the relative
energies (NEB, DFT, [synMTP], [eqMTP]), then the absolute energies in the same order.
An MTP column pair appears only if its option is given.
stdout: forward barrier (max - first) and reaction energy (last - first) of every curve.

neb_DFT.xyz: extxyz with the same ASE Atoms objects the DFT energies come from (one
  per OUTCAR, same order). It stores TOTEN as ``free_energy`` and energy(sigma->0) as
  ``energy``; forces and stress are included.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np

import matplotlib
from ase import Atoms
from ase.geometry import find_mic
from ase.io import read, write

from thermo_SiO2.io import read_sil

# Enable use on a compute node with no display server.
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# DEFAULT_NEB_DIR = Path(__file__).resolve().parent.parent / "d.neb"
DEFAULT_NEB_DIR = Path(__file__).resolve().parent / '../d.neb/'
DEFAULT_DFT_DIR = Path(__file__).resolve().parent / './calc'
OUTPUT_NAME = "fig_neb.pdf"
DFT_XYZ_NAME = "neb_DFT.xyz"
# E0_PATTERN = re.compile(r"energy\(sigma->0\)\s*=\s*([-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?)")


def _last_neb_profile(log_file: Path) -> tuple[list[float], list[float]]:
    """Read RDn/PEn values from the final row of the NEB master log."""
    header: list[str] | None = None
    lines = log_file.read_text(errors="replace").splitlines()

    for line in lines:
        fields = line.split()
        if "RD1" in fields and "PE1" in fields:
            header = fields

    if header is None:
        raise ValueError(f"No NEB RDn/PEn table found in {log_file}")

    fields = lines[-1].split() if lines else []
    if len(fields) != len(header):
        raise ValueError(
            f"Final line in {log_file} has {len(fields)} fields; "
            f"expected {len(header)} NEB fields"
        )

    try:
        values = dict(zip(header, map(float, fields)))
    except ValueError as exc:
        raise ValueError(f"Final line in {log_file} is not a numeric NEB row") from exc

    distances: list[float] = []
    energies: list[float] = []
    image = 1
    while f"RD{image}" in values and f"PE{image}" in values:
        distances.append(values[f"RD{image}"])
        energies.append(values[f"PE{image}"])
        image += 1

    if not distances:
        raise ValueError(f"Final line in {log_file} contains no RDn/PEn values")

    print(f'header = {header}')
    print(f'fields = {fields}')
    print(f'distances = {distances}')
    print(f'energies = {energies}')
    return distances, energies


def read_profile(neb_lmp_dir: Path) -> tuple[list[float], list[float]]:
    """Read the NEB master profile from the available ``log.lammps*`` files."""
    log_files = sorted(neb_lmp_dir.glob("log.lammps*"))
    if not log_files:
        raise FileNotFoundError(f"No log.lammps* files found in {neb_lmp_dir}")

    errors: list[str] = []
    for log_file in log_files:
        try:
            return _last_neb_profile(log_file)
        except ValueError as exc:
            errors.append(str(exc))
    raise ValueError("No usable NEB master log found:\n" + "\n".join(errors))


def read_grace_xyz_profile(xyz_file: Path) -> tuple[list[float], list[float]]:
    """Read normalized reaction coordinates and free energies from an NEB extxyz file."""
    images = read(xyz_file, index=":")
    if len(images) < 2:
        raise ValueError(f"{xyz_file} has {len(images)} image(s); need at least 2")

    # Cumulative 3N-dimensional distances with MIC, as in run_neb.print_neb_status.
    rd = [0.0]
    for previous, current in zip(images[:-1], images[1:]):
        displacement, _ = find_mic(
            current.positions - previous.positions, previous.cell, previous.pbc
        )
        rd.append(rd[-1] + float(np.linalg.norm(displacement)))
    if rd[-1] > 0:
        rd = [value / rd[-1] for value in rd]

    # extxyz free_energy= is stored on a SinglePointCalculator; force_consistent returns it.
    energies = [float(atoms.get_potential_energy(force_consistent=True)) for atoms in images]

    print(f'distances = {rd}')
    print(f'energies = {energies}')
    return rd, energies


def read_mtp_single_point_energies(cfg_file: Path, atom_symbols: str) -> list[float]:
    """Read MTP single-point energies from an MLIP cfg file."""
    cfgs = read_sil(cfg_file, atom_symbols=atom_symbols)
    energies = [float(cfg.get_potential_energy()) for cfg in cfgs]
    print(f'mtp_energies = {energies}')
    return energies


def read_dft_images(dft_dir: Path) -> list[Atoms]:
    """Read the last ionic step of each sorted ``calc/*/OUTCAR`` as ASE Atoms."""
    outcars = sorted(dft_dir.glob("*/OUTCAR"))
    if not outcars:
        raise FileNotFoundError(f"No calc/*/OUTCAR files found in {dft_dir}")
    return [read(outcar, index=-1, format="vasp-out") for outcar in outcars]


def plot_profile(
    distances: list[float], neb_energies: list[float], dft_energies: list[float], output: Path,
    neb_label: str = "LAMMPS NEB (r2SCAN-D4 syn-MTP)",
    syn_mtp_energies: list[float] | None = None,
    eq_mtp_energies: list[float] | None = None,
) -> None:
    """Write relative NEB, DFT, and optionally syn-/eq-MTP single-point energy profiles."""
    # (label, absolute energies, line style, color, column name)
    curves = [
        (neb_label, neb_energies, "o-", "C0", "NEB"),
        ("DFT single-point (r2SCAN-D4)", dft_energies, "s--", "C1", "DFT"),
    ]
    if syn_mtp_energies is not None:
        curves.append(("syn-MTP single-point (r2SCAN-D4)", syn_mtp_energies, "^:", "C2", "synMTP"))
    if eq_mtp_energies is not None:
        curves.append(("eq-MTP single-point (r2SCAN-D4)", eq_mtp_energies, "v-.", "C3", "eqMTP"))
    for label, energies, *_ in curves:
        if len(energies) != len(distances):
            raise ValueError(f"NEB has {len(distances)} images but {label} has {len(energies)}")

    relative = [np.asarray(energies) - energies[0] for _, energies, *_ in curves]

    fig, ax = plt.subplots(figsize=(5.0, 3.5), constrained_layout=True)
    for (label, _, style, color, _), rel in zip(curves, relative):
        ax.plot(distances, rel, style, color=color, linewidth=1.4, label=label)
    ax.axhline(0.0, color="0.6", linewidth=0.8, zorder=0)
    ax.set_xlabel("Reaction coordinate (normalized)")
    ax.set_ylabel("Energy (eV)")
    ax.margins(x=0.04, y=0.12)
    ax.tick_params(direction="in", top=True, right=True)
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)

    names = [curve[4] for curve in curves]
    header = "Reaction_coordinate_normalized  " + "  ".join(
        [f"{name}_rel_eV" for name in names] + [f"{name}_abs_eV" for name in names])
    columns = [distances, *relative, *[curve[1] for curve in curves]]
    np.savetxt('plot_data.txt', np.column_stack(columns), header=header, fmt='%.5f')

    for (label, *_), rel in zip(curves, relative):
        print(f"{label:35s} barrier {rel.max():.6f} eV  reaction {rel[-1]:.6f} eV")


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot LAMMPS NEB and DFT energy profiles.")
    neb_source = parser.add_mutually_exclusive_group()
    neb_source.add_argument("--neb-lmp-dir", type=Path, default=DEFAULT_NEB_DIR,
                            help="directory containing log.lammps* files (default: %(default)s)")
    neb_source.add_argument("--neb-grace-xyz", type=Path, default=None,
                            help="extxyz file with all GRACE NEB images; uses free_energy")
    parser.add_argument("--dft-dir", type=Path, default=DEFAULT_DFT_DIR,
                        help="directory containing ./*/OUTCAR files (default: %(default)s)")
    parser.add_argument("--single-point-syn-mtp", type=Path, default=None,
                        help="MLIP cfg file with syn-MTP single-point energies of the NEB images")
    parser.add_argument("--single-point-eq-mtp", type=Path, default=None,
                        help="MLIP cfg file with eq-MTP single-point energies of the NEB images")
    parser.add_argument("--atom-symbols", type=str, default="Si O H Al",
                        help="atom symbols passed to read_sil (default: %(default)s)")
    parser.add_argument("--output", type=Path, default=Path(OUTPUT_NAME),
                        help="output PDF path (default: %(default)s)")
    args = parser.parse_args()

    if args.neb_grace_xyz is not None:
        distances, energies = read_grace_xyz_profile(args.neb_grace_xyz)
        neb_label = "GRACE NEB"
    else:
        distances, energies = read_profile(args.neb_lmp_dir)
        neb_label = "LAMMPS NEB (r2SCAN-D4 syn-MTP)"
    dft_images = read_dft_images(args.dft_dir)
    # TOTEN (free energy), as in the mlip cfg.
    dft_energies = [float(atoms.get_potential_energy(force_consistent=True)) for atoms in dft_images]
    print(f'{dft_energies =}')
    write(DFT_XYZ_NAME, dft_images, format="extxyz")
    print(f"Wrote {DFT_XYZ_NAME} ({len(dft_images)} images)")
    syn_mtp_energies = None
    if args.single_point_syn_mtp is not None:
        syn_mtp_energies = read_mtp_single_point_energies(args.single_point_syn_mtp, args.atom_symbols)
    eq_mtp_energies = None
    if args.single_point_eq_mtp is not None:
        eq_mtp_energies = read_mtp_single_point_energies(args.single_point_eq_mtp, args.atom_symbols)
    plot_profile(distances, energies, dft_energies, args.output, neb_label=neb_label,
                 syn_mtp_energies=syn_mtp_energies, eq_mtp_energies=eq_mtp_energies)
    print(f"Wrote {args.output} ({len(energies)} images)")


if __name__ == "__main__":
    main()
