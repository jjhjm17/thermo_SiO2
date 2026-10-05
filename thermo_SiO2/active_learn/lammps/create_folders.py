#!/usr/bin/env python3
"""This script creates folders and prepares jobs for the active learning, for
lammps.

option 'variable_file', 'blanks', 'variables'
 The option should be turned on only when the keyword are present. So here r(param, 'variable_file',
  'in.file')
  If there is no 'blanks' variable, then and no 'variables' variable, the code should work without
  replacing blanks, without making errors.

When the code is changed, previous input files, which are not tested, should work, as much as possible.
Do not remove simply because it is not used in the current test files.

Keep this documentaion as it is hand-written.
"""
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
from thermo_SiO2.io import read_sil

def create_folders():
    """This function creates folders and prepares jobs."""

    if os.path.isfile('jobList'):
        os.remove('jobList')

    if use_initial_config:
        # structure_file = './hole_datafs/lammps.dataf_0'
        if type(initial_config) is str:
            structure_file = initial_config
        else:  # list
            structure_file = initial_config[0]
        # my_atoms = read(structure_file, format='lammps-data', style='atomic')
        # if symbols == 'Si O H Al':
        #     atom_symbol_tuple = Si_O_H_Al_atom_symbol_tuple_lammps
        # elif symbols == 'Si O Al':
        #     atom_symbol_tuple = Si_O_Al_atom_symbol_tuple_lammps
        # else:
        #     print("Error: the 'symbols' variable is not properly set.")
        #     sys.exit()
        # my_atoms = set_actual_atom_symbols(my_atoms,
        #         atom_symbol_tuple)
        my_atoms = read_sil(structure_file, atom_symbols=param.symbols)
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

        if hasattr(param, 'variable_file'):
            variable_file = param.variable_file
        elif os.path.isfile('in.file'):
            variable_file = 'in.file'
        else:
            variable_file = None

        if variable_file is not None:
            blanks = []
            variables = []

            if hasattr(param, 'blanks') != hasattr(param, 'variables'):
                raise ValueError("'blanks' and 'variables' must be defined together")
            if hasattr(param, 'blanks'):
                blanks.extend(param.blanks)
                variables.extend(str(value) for value in param.variables)

            if different_seeds:
                seed = get_lammps_random_seed(rng)
                blanks.append('xxxSEEDxxx')
                variables.append(str(seed))
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

            if blanks or variables:
                fill_blanks(file=variable_file, blanks=blanks,
                            variables=variables)

        # shutil.copy(f'../../hole_datafs/lammps.dataf_{index_seed}',
        #         'lammps.dataf')
        if use_initial_config:
            if type(initial_config) is str:
                shutil.copy(f'../../{initial_config}',
                        'lammps.dataf')
            else:  # list
                path_in = initial_config[index_seed]
                if path_in[0] == '.':  # relative path
                    path = f'../../{path_in}'
                else:
                    path = path_in
                shutil.copy(path,
                        'lammps.dataf')
        with open('../../jobList', 'a') as file:
            file.write(f'{os.getcwd()}\n')
        print(f'{calc_subfolder} ', end='')
        os.chdir('..')

    os.chdir('..')
    print()

if __name__ == "__main__":
    create_folders()
