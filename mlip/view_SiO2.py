"""This module has functions for viewing SiO2 configuration files, after 
changing the atom index to correct types of Si, O, ..."""

import argparse
from ase.visualize import view
from ase.io import read
# from .view_SiO2_mlip  import view_SiO2_mlip
# from .view_SiO2_lammps import view_SiO2_lammps
# from .view_SiO2_lammps_dump import view_SiO2_lammps_dump
from .read_mlip_cfg import read_mlip_cfg_mlippy
from .read_config import read_SiO2_dump
from ..util.SiO2_parameter import (Si_O_H_Al_atom_symbol_tuple_mlip,
        Si_O_H_Al_atom_symbol_tuple_lammps)
from ..util.util import set_actual_atom_symbols


def view_SiO2():
    """This function shows lammps or mlip SiO2 file."""
    # split = file.split('.')
    # if len(split) == 0:
    # suffix = file.split('.')[-1]
    # suffix_2nd = file.split('.')[-2]
    # if suffix == 'cfg' or suffix_2nd == 'cfg':

    parser = argparse.ArgumentParser(
        description='Shows SiO2 structures in lammps or mlip format.')
    parser.add_argument('--verbose', '-v', action='store_true')  # on/off
    parser.add_argument('--SiOAl', '-Al', action='store_true',
                       help='The order of atoms is Si O Al')
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
        # configs = read(file, format='lammps-dump-text', index=':')
        # for config in configs:
        #     config = set_actual_atom_symbols(config,
        #             Si_O_H_Al_atom_symbol_tuple_lammps)
        configs = read_SiO2_dump(file, index=':')
    else:  # lammps-data
        config = read(file, format='lammps-data', style='atomic')
        config = set_actual_atom_symbols(config,
                Si_O_H_Al_atom_symbol_tuple_lammps)
        configs = [config]


    if type(configs).__name__ == 'Atoms':
        configs = [configs]
        # If there is one configuration, convert to a list of configs.

    if args.SiOAl:
        #                                Si       O      H to Al
        atomic_number_to_symbol_SiOAl = {14:'Si', 8:'O', 1:'Al'}
        configs = [set_actual_atom_symbols(config,
                atomic_number_to_symbol_SiOAl) for config in configs]

    print(f'{len(configs)} configuration(s).')


    if args.verbose:
        for config in configs:
            print(config)
    view(configs)

if __name__ == '__main__':
    view_SiO2()
