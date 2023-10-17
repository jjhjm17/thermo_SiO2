#!/usr/bin/env python3

import sys
from ase.io import read
from ase.visualize import view
from ..util.util import set_actual_atom_symbols
from ..util.SiO2_parameter import Si_O_H_Al_atom_symbol_tuple_lammps

def view_SiO2_lammps(file):
    """This function shows lammps SiO2 file."""
    atoms_list = read(file, format='lammps-data', style='atomic')
    atoms_list = set_actual_atom_symbols(atoms_list,
            Si_O_H_Al_atom_symbol_tuple_lammps)
    for atoms in atoms_list:
        print(atoms)
    print(atoms_list)
    view(atoms_list)

if __name__ == '__main__':
    view_SiO2_lammps(sys.argv[1])
