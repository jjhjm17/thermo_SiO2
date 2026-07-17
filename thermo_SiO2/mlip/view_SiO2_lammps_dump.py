#!/usr/bin/env python3
"""This file reads a lammps dump file and shows SiO2."""

import sys
import ase.atoms
from ase.io import read
from ase.visualize import view
from ..util.util import set_actual_atom_symbols
from ..util.SiO2_parameter import Si_O_H_Al_atom_symbol_tuple_lammps

def view_SiO2_lammps_dump(file):
    """This function reads a lammps dump file and shows SiO2."""
    atoms_list = read(file, format='lammps-dump-text', index=':')
    # print(atoms_list)
    if type(atoms_list).__name__ == 'Atoms':
        atoms_list = [atoms_list]
        # Convert to a list of configs.
    for config in atoms_list:
        config = set_actual_atom_symbols(config,
                Si_O_H_Al_atom_symbol_tuple_lammps)
    for config in atoms_list:
        print(config)
    # print(atoms_list)
    view(atoms_list)

if __name__ == '__main__':
    view_SiO2_lammps(sys.argv[1])
