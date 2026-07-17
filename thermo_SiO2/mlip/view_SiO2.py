"""This module has functions for viewing SiO2 configuration files, after 
changing the atom index to correct types of Si, O, ..."""

import sys
import argparse
from ase.visualize import view
from ase.io import read
from .read import read_SiO2


def view_SiO2():
    """This function shows lammps or mlip SiO2 file.
    The default order is Si O H Al."""
    # split = file.split('.')
    # if len(split) == 0:
    # suffix = file.split('.')[-1]
    # suffix_2nd = file.split('.')[-2]
    # if suffix == 'cfg' or suffix_2nd == 'cfg':

    parser = argparse.ArgumentParser(
        description='Shows SiO2 structures in lammps or mlip format.')
    parser.add_argument("files", nargs="+", help="One or more structure files to view.")
    parser.add_argument('--verbose', '-v', action='store_true')  # on/off
    parser.add_argument('--SiOAl', '-Al', '-al', action='store_true',
                       help='The order of atoms is Si O Al')
    parser.add_argument('--SiOY', '-Y', '-y', action='store_true',
                       help='The order of atoms is Si O Y')
    parser.add_argument('--wrap', '-w', action='store_true',
                       help='Wrap atoms')
    parser.add_argument('--atom_symbols', '-s', type=str,
                       help="atom_symbols, ex) 'O H'", default='Si O H Al')
    args = parser.parse_args()
    files = args.files
    if args.SiOAl:
        atom_symbols = 'Si O Al'
    elif args.SiOY:
        atom_symbols = 'Si O Y'
    else:
        # atom_symbols = 'Si O H Al'
        atom_symbols = args.atom_symbols

    all_configs = []
    silent = True
    for ind, file in enumerate(files):
        configs = read_SiO2(file, atom_symbols, silent=silent)
        # print(f'{configs = }')
        all_configs.extend(configs)
    print(f'{len(all_configs)} snapshot(s).')
    print(f'{all_configs[0] = }')

    if args.verbose:
        for i_cfg, config in enumerate(all_configs):
            print(f'{i_cfg}  {config}')
    if args.wrap:
        for config in all_configs:
            config = config.wrap()
    view(all_configs)

if __name__ == '__main__':
    view_SiO2()
