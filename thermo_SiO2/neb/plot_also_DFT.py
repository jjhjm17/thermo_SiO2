"""
This script reads the final energies of (for example) ../d.neb/log.lammps.* files and plot
Energy (eV) vs Reaction coordinate (Ang) to file fig_neb.pdf.

Also plot DFT energies in (for example) ../f.DFT/calc/*/OUTCAR folders, sorted. Read energy(sigma->0) (E0).

How to find the last line?
After getting the header, check the last line and see if it is valid, that is, it has the same number of whitespace-
  separated fields.

Keep comments in this comment block, these are written by a person.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import numpy as np

import matplotlib

# Enable use on a compute node with no display server.
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# DEFAULT_NEB_DIR = Path(__file__).resolve().parent.parent / "d.neb"
DEFAULT_NEB_DIR = Path(__file__).resolve().parent / '../d.neb/'
DEFAULT_DFT_DIR = Path(__file__).resolve().parent / './calc'
OUTPUT_NAME = "fig_neb.pdf"
# E0_PATTERN = re.compile(r"energy\(sigma->0\)\s*=\s*([-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?)")
E0_PATTERN = re.compile(r"energy\(sigma->0\) =\s*(\S*)")
# \s: empty space,  \S: not empty space, *: zero or many


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


def read_profile(neb_dir: Path) -> tuple[list[float], list[float]]:
    """Read the NEB master profile from the available ``log.lammps*`` files."""
    log_files = sorted(neb_dir.glob("log.lammps*"))
    if not log_files:
        raise FileNotFoundError(f"No log.lammps* files found in {neb_dir}")

    errors: list[str] = []
    for log_file in log_files:
        try:
            return _last_neb_profile(log_file)
        except ValueError as exc:
            errors.append(str(exc))
    raise ValueError("No usable NEB master log found:\n" + "\n".join(errors))


def _last_dft_energy(outcar: Path) -> float:
    """Return the final VASP energy(sigma->0) value from an OUTCAR."""
    matches = E0_PATTERN.findall(outcar.read_text(errors="replace"))
    if not matches:
        raise ValueError(f"No energy(sigma->0) value found in {outcar}")
    return float(matches[-1])


def read_dft_energies(dft_dir: Path) -> list[float]:
    """Read final E0 values from sorted ``calc/*/OUTCAR`` calculations."""
    outcars = sorted(dft_dir.glob("*/OUTCAR"))
    if not outcars:
        raise FileNotFoundError(f"No calc/*/OUTCAR files found in {dft_dir}")
    return [_last_dft_energy(outcar) for outcar in outcars]


def plot_profile(
    distances: list[float], neb_energies: list[float], dft_energies: list[float], output: Path
) -> None:
    """Write relative LAMMPS NEB and DFT energy profiles."""
    if len(distances) != len(dft_energies):
        raise ValueError(
            f"NEB has {len(distances)} images but DFT has {len(dft_energies)} OUTCAR files"
        )

    relative_neb_energies = [energy - neb_energies[0] for energy in neb_energies]
    relative_dft_energies = [energy - dft_energies[0] for energy in dft_energies]
    fig, ax = plt.subplots(figsize=(5.0, 3.5), constrained_layout=True)
    ax.plot(distances, relative_neb_energies, "o-", color="C0", linewidth=1.6, label="LAMMPS NEB (r2SCAN-D4 syn-MTP)")
    ax.plot(distances, relative_dft_energies, "s--", color="C1", linewidth=1.4, label="DFT single-point (r2SCAN-D4)")
    ax.axhline(0.0, color="0.6", linewidth=0.8, zorder=0)
    ax.set_xlabel("Reaction coordinate (normalized)")
    ax.set_ylabel("Energy (eV)")
    ax.margins(x=0.04, y=0.12)
    ax.tick_params(direction="in", top=True, right=True)
    ax.legend(frameon=False)
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)

    np.savetxt('plot_data.txt',
               np.column_stack([distances, relative_neb_energies, relative_dft_energies]),
               header='Reaction coord. (normalized), relative NEB E (eV), DFT E (eV)',
               fmt='%.6g')


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot LAMMPS NEB and DFT energy profiles.")
    parser.add_argument("--neb-dir", type=Path, default=DEFAULT_NEB_DIR,
                        help="directory containing log.lammps* files (default: %(default)s)")
    parser.add_argument("--dft-dir", type=Path, default=DEFAULT_DFT_DIR,
                        help="directory containing ./*/OUTCAR files (default: %(default)s)")
    parser.add_argument("--output", type=Path, default=Path(OUTPUT_NAME),
                        help="output PDF path (default: %(default)s)")
    args = parser.parse_args()

    distances, energies = read_profile(args.neb_dir)
    dft_energies = read_dft_energies(args.dft_dir)
    print(f'{dft_energies =}')
    plot_profile(distances, energies, dft_energies, args.output)
    print(f"Wrote {args.output} ({len(energies)} images)")
    print(f"Forward barrier: {max(energies) - energies[0]:.6f} eV")
    print(f"Reaction energy: {energies[-1] - energies[0]:.6f} eV")

    np.savetxt('plot_data.txt', np.column_stack([distances, energies, dft_energies]),
               header='Reaction_coordiate_normalized  MTP_NEB_eV DFT_eV')


if __name__ == "__main__":
    main()
