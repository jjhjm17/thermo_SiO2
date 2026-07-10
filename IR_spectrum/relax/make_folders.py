#!/usr/bin/env python3
"""This script makes folders and prepares jobs for the melt quench process."""

import os
import sys
from shutil import copy
import numpy as np
from ase.visualize import view
from ase.io import read, write
from ...util.util import get_lammps_random_seed, set_actual_atom_symbols
from ...util.SiO2_parameter import (Si_O_Al__Z_of_type_lammps,
                                    Si_O_H_Al__Z_of_type_lammps,
                                    Si_O_Al_atom_symbol_tuple_lammps,
                                    Si_O_H_Al_atom_symbol_tuple_lammps,
                                    spec_order_POSCAR_to_LAMMPS)
from a_parameters import calc_folder, random_seed, structures, symbols
import a_parameters
import a_parameters as param


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
        if hasattr(a_parameters, 'format'):
            format = a_parameters.format
        else:
            format == 'lammps-data'  # default

        if symbols == 'Si O H Al':
            atom_symbol_tuple = Si_O_H_Al_atom_symbol_tuple_lammps
            Z_of_type = Si_O_H_Al__Z_of_type_lammps
        elif symbols == 'Si O Al':
            atom_symbol_tuple = Si_O_Al_atom_symbol_tuple_lammps
            Z_of_type = Si_O_Al__Z_of_type_lammps
        else:
            print("Error: the 'symbols' variable is not properly set.")
            sys.exit()

        if config_file.startswith('/'):
            cfg_path = config_file
        else:
            cfg_path = f'../{config_file}'

        if format == 'lammps-dump-text':
            # my_atoms = read_SiO2_dump(f'../{config_file}', index=-1)
            my_atoms = read(cfg_path, format=format)
            my_atoms = set_actual_atom_symbols(my_atoms,
                atom_symbol_tuple)
        elif format == 'lammps-data':
            if hasattr(param, 'atom_style_in'):
                atom_style = param.atom_style_in
            else:
                atom_style='atomic'
            # lammps-data
            my_atoms = read(cfg_path, format=format, atom_style=atom_style,
                            Z_of_type=Z_of_type)
        my_atoms.wrap()
        print(f'  {my_atoms}')
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
        # os.symlink('../../mlip.ini', 'mlip.ini')
        # copy(f'../../{config_file}', 'input_structure.dataf')
        # write('input_structure.dataf', my_atoms, format='lammps-data',
        #       spec_order=['Si', 'O', 'H', 'Al'])
        if hasattr(param, 'charge') and param.charge:
            for atom in my_atoms:
                atom.charge = 0
            write('coords.dataf', my_atoms, format='lammps-data',
                  specorder=spec_order_POSCAR_to_LAMMPS, atom_style='charge')
        else:
            write('coords.dataf', my_atoms, format='lammps-data',
                  specorder=spec_order_POSCAR_to_LAMMPS)
        with open('../../jobList', 'a') as file:
            file.write(f'{os.getcwd()}\n')
        os.chdir('..')
    os.chdir('..')


if __name__ == "__main__":
    make_folders()
