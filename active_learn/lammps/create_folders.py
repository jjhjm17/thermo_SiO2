#!/usr/bin/env python3
"""This script creates folders and prepares jobs for the active learning, for
lammps."""
import os
import sys
import shutil
import numpy as np
from ase.visualize import view
from ase.io import read
from a_parameters import calc_folder, num_seeds, use_initial_config, initial_config, symbols, different_seeds, template_folder
import a_parameters as param
from ...util.util import (get_lammps_random_seed, set_actual_atom_symbols,
                          fill_blanks)
from ...util.SiO2_parameter import (Si_O_H_Al_atom_symbol_tuple_lammps,
                                    Si_O_Al_atom_symbol_tuple_lammps)

def create_folders():
    """This function creates folders and prepares jobs."""

    if os.path.isfile('jobList'):
        os.remove('jobList')

    if use_initial_config:
        # structure_file = './hole_datafs/lammps.dataf_0'
        structure_file = initial_config
        my_atoms = read(structure_file, format='lammps-data', style='atomic')
        if symbols == 'Si O H Al':
            atom_symbol_tuple = Si_O_H_Al_atom_symbol_tuple_lammps
        elif symbols == 'Si O Al':
            atom_symbol_tuple = Si_O_Al_atom_symbol_tuple_lammps
        else:
            print("Error: the 'symbols' variable is not properly set.")
            sys.exit()
        my_atoms = set_actual_atom_symbols(my_atoms,
                atom_symbol_tuple)
        print(my_atoms)
        view(my_atoms)
        print('Please check the structure visually.')

    template_folder_abs = os.path.abspath(template_folder)
    if not os.path.exists(calc_folder):
        os.mkdir(calc_folder)
    os.chdir(calc_folder)

    print(f'folders made: {calc_folder}/ ', end='')

    if different_seeds:
        if hasattr(param, 'seed_of_seeds'):
            rng = np.random.default_rng(param.seed_of_seeds)
        else:
            rng = None

    for index_seed in range(num_seeds):
        calc_subfolder = f'{str(index_seed).zfill(3)}'
        # shutil.copytree(f'../template/{template_folder}', calc_subfolder, symlinks=True)
        shutil.copytree(template_folder_abs, calc_subfolder, symlinks=True)
        os.chdir(calc_subfolder)

        if different_seeds:
            if hasattr(param, 'variable_file'):
                variable_file = param.variable_file
            else:
                variable_file = 'in.file'
            with open(variable_file, 'r') as fin:
                in_file_lines = fin.readlines()
            seed = get_lammps_random_seed(rng)
            # with open('in.file_new', 'w') as out_file:
            #     for line in in_file_lines:
            #         line = line.replace('xxxSEEDxxx', f'{seed}')
            #         if hasattr(param, 'almtp'):
            #             line = line.replace('xxx__almpt__xxx',
            #                                 param.almtp)
            #         out_file.write(line)
            # os.remove('in.file')
            # os.rename('in.file_new', 'in.file')
            blanks = ['xxxSEEDxxx']
            variables = [str(seed)]
            if hasattr(param, 'almtp'):
                blanks.append('xxx__almtp__xxx')
                if param.almtp.startswith('/'):
                    almtp_path = param.almtp
                else:
                    almtp_path = f'../../{param.almtp}'
                variables.append(almtp_path)
            if (hasattr(param, 'fill_blank_index') and
                param.fill_blank_index):
                blanks.append('xxx__index__xxx')
                variables.append(f'{index_seed}')
            fill_blanks(file=variable_file, blanks=blanks, variables=variables)

        # shutil.copy(f'../../hole_datafs/lammps.dataf_{index_seed}',
        #         'lammps.dataf')
        if use_initial_config:
            shutil.copy(f'../../{initial_config}',
                    'lammps.dataf')
        with open('../../jobList', 'a') as file:
            file.write(f'{os.getcwd()}\n')
        print(f'{calc_subfolder} ', end='')
        os.chdir('..')

    os.chdir('..')
    print()

if __name__ == "__main__":
    create_folders()
