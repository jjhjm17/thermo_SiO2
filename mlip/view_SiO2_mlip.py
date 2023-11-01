#!/usr/bin/env python3

import sys
from .read_mlip_cfg import read_mlip_cfg_positions
from ase.visualize import view
from ..util.util import set_actual_atom_symbols
from ..util.SiO2_parameter import Si_O_H_Al_atom_symbol_tuple_mlip

def view_SiO2_mlip(file):
    """This function shows lammps SiO2 file."""
    atoms_list = read_mlip_cfg_positions(file,
            atom_symbol_tuple=Si_O_H_Al_atom_symbol_tuple_mlip)
    if type(atoms_list).__name__ == 'Atoms':
        atoms_list = [atoms_list]
        # Convert to a list of configs.
    for atoms in atoms_list:
        print(atoms)
    view(atoms_list)

if __name__ == '__main__':
    view_SiO2_mlip(sys.argv[1])
