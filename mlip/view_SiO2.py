#!/usr/bin/env python3

import sys
import argparse
from ase.visualize import view
# from .view_SiO2_mlip  import view_SiO2_mlip
# from .view_SiO2_lammps import view_SiO2_lammps
# from .view_SiO2_lammps_dump import view_SiO2_lammps_dump
from .read_mlip_cfg import read_mlip_cfg_positions, read_mlip_cfg_mlippy
from ..util.SiO2_parameter import (Si_O_H_Al_atom_symbol_tuple_mlip,
        Si_O_H_Al_atom_symbol_tuple_lammps)


def view_SiO2():
    """This function shows lammps or mlip SiO2 file."""
    # split = file.split('.')
    # if len(split) == 0:
    # suffix = file.split('.')[-1]
    # suffix_2nd = file.split('.')[-2]
    # if suffix == 'cfg' or suffix_2nd == 'cfg':

    parser = argparse.ArgumentParser(
        description='Shows SiO2 structures in lammps or mlip format.')
    parser.add_argument('--verbose', '-v', action='store_true')
    parser.add_argument("file")
    args = parser.parse_args()
    file = args.file

    if '.cfg' in file:
        # configs: configurations
        # configs = read_mlip_cfg_positions(file,
        #         atom_symbol_tuple=Si_O_H_Al_atom_symbol_tuple_mlip)
        configs = read_mlip_cfg_mlippy(file,
                atom_symbol_tuple=Si_O_H_Al_atom_symbol_tuple_mlip)
    elif 'dump' in file:
        configs = read(file, format='lammps-dump-text', index=':')
        for config in configs:
            config = set_actual_atom_symbols(config,
                    Si_O_H_Al_atom_symbol_tuple_lammps)
    else:  # lammps-data
        configs = read(file, format='lammps-data', style='atomic')
        configs = set_actual_atom_symbols(configs,
                Si_O_H_Al_atom_symbol_tuple_lammps)

    if type(configs).__name__ == 'Atoms':
        configs = [configs]
        # If there is one configuration, convert to a list of configs.
    print(f'{len(configs)} configurations.')
    if args.verbose:
        for config in configs:
            print(config)
    view(configs)

if __name__ == '__main__':
    view_SiO2()
