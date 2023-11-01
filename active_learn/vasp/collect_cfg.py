#!/usr/bin/env python3
"""This script collects calculated configuration to cfg file."""

import os
import sys
from parameters import mlip_cfg_file, calc_folder, output_file
from ...util.util import shell
from ...mlip.read_mlip_cfg import read_mlip_cfg
# from ...mlip.read_mlip_cfg import set_atom_symbol
from ...util.SiO2_parameter import Si_O_H_Al_atom_symbol_tuple_mlip

def collect_cfg():
    """This function collects calculated configurations."""

    if os.path.exists(output_file):
        print(f'Error: output file {output_file} exists.')
        sys.exit()

    # atoms = ase_loadcfgs(mlip_cfg_file)
    # set_atom_symbol(atoms, Si_O_H_Al_atom_symbol_tuple_mlip)
    atoms_and_forces = read_mlip_cfg(mlip_cfg_file,
            atom_symbol_tuple=Si_O_H_Al_atom_symbol_tuple_mlip,
            sort_method='POTCAR_order')

    os.chdir(calc_folder)

    for i_structure in range( len( atoms_and_forces)):
        # i_; index
        # struct_str = '{:03d}'.format(i_structure)    # string, padded with 0
        struct_str = f'{i_structure:03d}'    # string, padded with 0
        folder = struct_str
        print(f'directory = {folder}   ', end='')

        os.chdir(folder)
        shell(f'mlp3ser convert OUTCAR ../../{output_file} --append '
            '--input_format=outcar')
        os.chdir('..')
    os.chdir('..')

    atoms_and_forces_out = read_mlip_cfg(output_file,
            atom_symbol_tuple=Si_O_H_Al_atom_symbol_tuple_mlip,
            sort_method='POTCAR_order')
    print(f'Total {len(atoms_and_forces_out)} configurations are collected in '
          f'this file: {output_file}')

if __name__ == '__main__':
    collect_cfg()
