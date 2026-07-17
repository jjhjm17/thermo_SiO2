#!/usr/bin/env python3
"""This script creates folders and prepares jobs for the melt quench process."""

import os
import sys
import numpy as np
from ase.visualize import view
from ase.io import read
from ...util.util import get_lammps_random_seed
from a_parameters import num_seeds, calc_folder, input_file


def create_folders():
    """This function creates folders and prepares jobs for melt quench and
    equilibration."""

    if os.path.isfile('jobList'):
        os.remove('jobList')

    structure_file = 'cristo_192.dataf'
    my_atoms = read(structure_file, format='lammps-data', style='atomic')
    # num_atoms = len(my_atoms)
    # my_atoms.set_chemical_symbols(f'{atom_symbol}{num_atoms}')
    print('Please check the structure visually.')
    view(my_atoms)

    if not os.path.exists(calc_folder):
        os.mkdir(calc_folder)
    os.chdir(calc_folder)
    for index_seed in range(1, num_seeds+1):  # 1 to num_seeds
        calc_subfolder = f'seed_{str(index_seed).zfill(3)}'
        print(calc_subfolder)
        os.mkdir(calc_subfolder)
        os.chdir(calc_subfolder)
        seed = get_lammps_random_seed()
        with open('../../in.file', 'r') as in_file:
        # with open(f'../../{input_file}', 'r') as in_file:
            in_file_lines = in_file.readlines()

        with open('./in.file', 'w') as out_file:
            for line in in_file_lines:
                line = line.replace('xxxSEEDxxx', f'{seed}')
                out_file.write(line)
        os.symlink('../../SiO2_amorphous_lda_level22.mtp',
                'SiO2_amorphous_lda_level22.mtp')
        os.symlink('../../mlip.ini', 'mlip.ini')
        os.symlink(f'../../{structure_file}', structure_file)
        with open('../../jobList', 'a') as file:
            file.write(f'{os.getcwd()}\n')
        os.chdir('..')
    os.chdir('..')


if __name__ == "__main__":
    create_folders()
