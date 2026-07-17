"""This module contains constants and methods for SiO2 parameters.
"""
from ase.data import chemical_symbols

POTCAR_setup = {'Si': 'Si', 'O': 'O', 'H': 'H', 'Al': 'Al', 'Y': 'Y_sv'}

# The variables for orders will be replaced by the class Atom_order below.

# Si O H Al
Si_O_H_Al_atom_symbol_tuple_lammps = {1:'Si', 2: 'O', 3: 'H', 4: 'Al'}
# from POTCAR

Si_O_H_Al_atom_symbol_tuple_mlip = {0:'Si', 1: 'O', 2: 'H', 3: 'Al'}

atomic_symbol_to_number = {"H": 1, "O": 8, "Si": 14, "Al": 13}

                              #  Si     O     H     Al
atomic_number_to_POTCAR_order = {14: 0, 8: 1, 1: 2, 13: 3}

spec_order_POSCAR_to_LAMMPS = ['Si', 'O', 'H', 'Al']

#                              Si    O    H    Al
Si_O_H_Al__Z_of_type_lammps = {1:14, 2:8, 3:1, 4:13}


# Si O Al
Si_O_Al_atom_symbol_tuple_mlip = {0:'Si', 1: 'O', 2: 'Al'}
Si_O_Al_atom_symbol_tuple_lammps = {1:'Si', 2: 'O', 3: 'Al'}

#                              Si    O    Al
Si_O_Al__Z_of_type_lammps = {1:14, 2:8, 3:13}
                                       #  Si     O     Al
Si_O_Al_atomic_number_to_POTCAR_order = {14: 0, 8: 1, 13: 2}


class Atom_order:
    """ symbols: atomic symbols, ex) 'Si O H Al'

    symbols_list: ex) ['Si', 'O', 'H', 'Al']
    indices; ex) [0, 1, 2, 3]
    Z_values: atomic numbers starting from 1, ex) [14, 8, 1, 13]

    """
    def __init__(self, symbols):
        self.symbol_list = symbols.split()
        self.indices = list(range(len(self.symbol_list)))
        self.indices_1 = list(range(1, len(self.symbol_list)+1))

        # Z_values = []
        # for symbol in self.symbol_list:
        #     Z_value = chemical_symbols.index(symbol)
        #     Z_values.append(Z_value)
        # self.Z_values = Z_values
        self.Z_values = [chemical_symbols.index(symbol) for symbol in
                         self.symbol_list]

    def atomic_number_to_POTCAR_order(self):
        """ Ex) #  Si     O     H     Al
                  {14: 0, 8: 1, 1: 2, 13: 3}"""
        return dict(zip(self.Z_values, self.indices))

    def atom_symbol_tuple_mlip(self):
        """
        atom_symbol_tuple_mlip:
            ex) {0:'Si', 1: 'O', 2: 'H', 3: 'Al'}"""
        return dict(zip(self.indices, self.symbol_list))

    def atom_symbol_tuple_lammps(self):
        """
        atom_symbol_tuple_lammps:
            ex) {1:'Si', 2: 'O', 3: 'H', 4: 'Al'}"""
        return dict(zip(self.indices_1, self.symbol_list))

    def Z_of_type_lammps(self):
        """                       Si    O    H    Al
            Z_of_type_lammps = {1:14, 2:8, 3:1, 4:13} """

        return dict(zip(self.indices_1, self.Z_values))
