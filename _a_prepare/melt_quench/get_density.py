#!/usr/bin/env python3
"""This script obtains densities of SiO2."""

import os
import numpy as np
from ase.io import read
from a_parameters import num_seeds, calc_folder
from ...util.util import read_output


def get_density():
    """This function cleans calculation files."""

    os.chdir(calc_folder)
    densities = []
    print('Folder   number of configs')
    for index_seed in range(1, num_seeds+1):  # 1 to num_seeds
        calc_subfolder = f'seed_{str(index_seed).zfill(3)}'
        os.chdir(calc_subfolder)
        num_configs = read_output('grep TIMESTEP dump_atom | wc --line')
        print(f'{calc_subfolder} {num_configs}')
        config = read('dump_atom', index='0', format='lammps-dump-text')
        for atom in config:
            if atom.symbol == 'H':
                atom.symbol = 'Si'
            elif atom.symbol == 'He':
                atom.symbol = 'O'
        amu_per_ang_3__to__g_per_cm_3 = 1.6605391
        # units u/angstrom^3 g/cm^3
        density = (np.sum(config.get_masses()) / config.get_volume() *
                amu_per_ang_3__to__g_per_cm_3)
        densities.append(density)
        os.chdir('..')
    os.chdir('..')
    avg_density = np.average(densities)
    std = np.std(densities, ddof=1)
    print(f'Density = {avg_density:.4f} +- {std:.4f} (std) g/cm^3 '
          f'(n={len(densities)})')

if __name__ == "__main__":
    get_density()
