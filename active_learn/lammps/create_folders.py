#!/usr/bin/env python3
"""This script creates folders and prepares jobs for the active learning, for
lammps."""
import os
import shutil
from ase.visualize import view
from ase.io import read
from ...util.util import get_lammps_random_seed, set_actual_atom_symbols
from ...util.SiO2_parameter import Si_O_H_Al_atom_symbol_tuple_lammps
from a_parameters import calc_folder, num_seeds

def create_folders():
    """This function creates folders and prepares jobs."""

    if os.path.isfile('jobList'):
        os.remove('jobList')

    structure_file = './hole_datafs/lammps.dataf_0'
    my_atoms = read(structure_file, format='lammps-data', style='atomic')
    my_atoms = set_actual_atom_symbols(my_atoms,
            Si_O_H_Al_atom_symbol_tuple_lammps)
    print(my_atoms)
    view(my_atoms)
    print('Please check the structure visually.')

    if not os.path.exists(calc_folder):
        os.mkdir(calc_folder)
    os.chdir(calc_folder)

    print(f'folders made: {calc_folder}/ ', end='')

    for index_seed in range(num_seeds):
        calc_subfolder = f'{str(index_seed).zfill(3)}'
        shutil.copytree('../template/0', calc_subfolder, symlinks=True)
        os.chdir(calc_subfolder)

        # with open('in.file', 'r') as fin:
        #     in_file_lines = fin.readlines()
        # seed = get_lammps_random_seed()
        # with open('in.file_new', 'w') as out_file:
        #     for line in in_file_lines:
        #         line = line.replace('xxxSEEDxxx', f'{seed}')
        #         out_file.write(line)
        # os.remove('in.file')
        # os.rename('in.file_new', 'in.file')

        shutil.copy(f'../../hole_datafs/lammps.dataf_{index_seed}',
                'lammps.dataf')
        with open('../../jobList', 'a') as file:
            file.write(f'{os.getcwd()}\n')
        print(f'{calc_subfolder} ', end='')
        os.chdir('..')

    os.chdir('..')
    print()

if __name__ == "__main__":
    create_folders()
