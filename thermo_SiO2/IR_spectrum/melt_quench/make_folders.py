#!/usr/bin/env python3
"""This script makes folders and prepares jobs for the melt quench process."""

import os
import sys
from shutil import copy
import numpy as np
from ase.visualize import view
from ase.io import read
from ...util.util import get_lammps_random_seed
from ...util.SiO2_parameter import Si_O_Al__Z_of_type_lammps
from a_parameters import calc_folder, random_seed, structures


def make_folders():
    """This function makes folders and prepares jobs for melt quench and
    equilibration."""

    if os.path.isfile('jobList'):
        os.remove('jobList')

    if not os.path.exists(calc_folder):
        os.mkdir(calc_folder)
    os.chdir(calc_folder)
    for i_config, config_file in enumerate(structures):
        print(f'index = {i_config}, {config_file}')
        my_atoms = read(f'../{config_file}', format='lammps-data',
                        style='atomic', Z_of_type=Si_O_Al__Z_of_type_lammps)
        print(f'  {my_atoms}')
        # num_atoms = len(my_atoms)
        # my_atoms.set_chemical_symbols(f'{atom_symbol}{num_atoms}')
        if i_config == 0:
            print('Please check the first structure visually.')
            # shell(f'view_SiO2 ../{config_file}')
            view(my_atoms)

        calc_subfolder = f'struct_{str(i_config).zfill(3)}'
        print(f'{calc_subfolder}\n')
        os.mkdir(calc_subfolder)
        os.chdir(calc_subfolder)
        if random_seed:
            seed = get_lammps_random_seed()
        with open('../../in.file', 'r') as in_file:
            in_file_lines = in_file.readlines()

        with open('./in.file', 'w') as out_file:
            for line in in_file_lines:
                if random_seed:
                    line = line.replace('xxxSEEDxxx', f'{seed}')
                out_file.write(line)
        os.symlink('../../mlip.ini', 'mlip.ini')
        copy(f'../../{config_file}', 'input_structure.dataf')
        with open('../../jobList', 'a') as file:
            file.write(f'{os.getcwd()}\n')
        os.chdir('..')
    os.chdir('..')


if __name__ == "__main__":
    make_folders()
