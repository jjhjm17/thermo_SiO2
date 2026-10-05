#!/usr/bin/env python3
"""Obtain densities of SiO2 configurations in g/cm^3."""

from pathlib import Path

import numpy as np
from ase.io import read


AMU_PER_ANGSTROM3_TO_G_PER_CM3 = 1.6605391


def get_density_cfg(cfg):
    """Return the density of an ASE Atoms configuration in g/cm^3.

    Atomic masses are in atomic mass units and the cell volume is in
    angstrom^3. The configuration must have its actual chemical species.
    """
    volume = cfg.get_volume()
    if not np.isfinite(volume) or volume <= 0:
        raise ValueError("Configuration must have a positive, finite cell volume")
    return float(np.sum(cfg.get_masses()) / volume *
                 AMU_PER_ANGSTROM3_TO_G_PER_CM3)


def get_density():
    """Report the mean and sample standard deviation across seed densities."""
    from a_parameters import num_seeds, calc_folder

    densities = []
    print('Folder   number of configs')
    for index_seed in range(1, num_seeds + 1):
        calc_subfolder = f'seed_{index_seed:03d}'
        dump_path = Path(calc_folder) / calc_subfolder / 'dump_atom'
        with dump_path.open() as dump_file:
            num_configs = sum(line.startswith('ITEM: TIMESTEP')
                              for line in dump_file)
        print(f'{calc_subfolder} {num_configs}')

        # Preserve the existing first-frame selection. LAMMPS types 1 and 2
        # represent Si and O, respectively.
        config = read(dump_path, index=0, format='lammps-dump-text',
                      specorder=['Si', 'O'])
        densities.append(get_density_cfg(config))

    if not densities:
        raise ValueError("No seeds available for density calculation")
    avg_density = np.average(densities)
    std = np.std(densities, ddof=1) if len(densities) > 1 else float('nan')
    print(f'Density = {avg_density:.4f} +- {std:.4f} (std) g/cm^3 '
          f'(n={len(densities)})')


if __name__ == "__main__":
    get_density()
