"""This script contains functions for reading SiO2 configurations."""

import sys
from ase.io import read
from .read_mlip_cfg import read_mlip_cfg_mlippy
from ..util.SiO2_parameter import Atom_order


def read_SiO2(file, atom_symbols, index=':'):
    """This function reads lammps or mlip SiO2 file generally."""

    order = Atom_order(atom_symbols)

    if '.cfg' in file:
        # configs: configurations
        atom_symbol_tuple = order.atom_symbol_tuple_mlip()
        configs = read_mlip_cfg_mlippy(file,
                atom_symbol_tuple=atom_symbol_tuple, index=index)
    elif 'dump' in file:
        specorder = atom_symbols.split()
        configs = read(file, index=index, format='lammps-dump-text', specorder=specorder)
    elif 'dataf' in file:  # lammps-data
        # config = read(file, format='lammps-data', atom_style='atomic')
        Z_of_type = order.Z_of_type_lammps()
        configs = read(file, index=index, format='lammps-data',
                       Z_of_type=Z_of_type, atom_style='atomic')
        # print(f'{config = }')
        # configs = [config]
    elif file.endswith('.db'):  # ase db
        print(f'{file = }')
        print(f'{index = }')
        configs = read(file, index=index)
    else:
        print('Error: unknown file format.')
        sys.exit()

    if type(configs).__name__ == 'Atoms':
        configs = [configs]
        # If there is one configuration, convert to a list of configs.

    print(f'{len(configs)} configuration(s).')
    print(f'{configs[0] = }')
    return configs
