"""
This script reads the final energies of (given folder, ex. ../d.neb) / log.lammps.* files and plot
Energy (eV) vs Reaction coordinate (Ang) to file fig_neb.pdf.

How to find the last line?
After getting the header, check the last line and see if it is valid, that is, it has the same number of whitespace-
  separated fields.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
import numpy as np

# Enable use on a compute node with no display server.
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# DEFAULT_NEB_DIR = Path(__file__).resolve().parent.parent / "d.neb"
DEFAULT_NEB_DIR = '.'
OUTPUT_NAME = "fig_neb.pdf"


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
    """Read the converged climbing-image NEB profile from ``log.lammps``."""
    return _last_neb_profile(neb_dir / "log.lammps")


def plot_profile(distances: list[float], energies: list[float], output: Path) -> None:
    """Write the NEB energy profile relative to the initial image."""
    relative_energies = [energy - energies[0] for energy in energies]
    fig, ax = plt.subplots(figsize=(5.0, 3.5), constrained_layout=True)
    ax.plot(distances, relative_energies, "o-", color="C0", linewidth=1.6)
    ax.axhline(0.0, color="0.6", linewidth=0.8, zorder=0)
    ax.set_xlabel("Reaction coordinate (normalized)")
    ax.set_ylabel("Energy (eV)")
    ax.margins(x=0.04, y=0.12)
    ax.tick_params(direction="in", top=True, right=True)
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot a converged LAMMPS NEB energy profile.")
    parser.add_argument("--neb-dir", type=Path, default=DEFAULT_NEB_DIR,
                        help="directory containing log.lammps (default: %(default)s)")
    parser.add_argument("--output", type=Path, default=Path(OUTPUT_NAME),
                        help="output PDF path (default: %(default)s)")
    args = parser.parse_args()

    distances, energies = read_profile(args.neb_dir)
    plot_profile(distances, energies, args.output)
    print(f"Wrote {args.output} ({len(energies)} images)")
    print(f"Forward barrier: {max(energies) - energies[0]:.6f} eV")
    print(f"Reaction energy: {energies[-1] - energies[0]:.6f} eV")

    np.savetxt('plot_data.txt', np.column_stack([distances, energies]),
               header='Reaction_coordiate_normalized  MTP_NEB_eV', fmt='%.5f')

if __name__ == "__main__":
    main()
