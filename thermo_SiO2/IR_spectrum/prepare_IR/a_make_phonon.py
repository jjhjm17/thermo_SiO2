#!/usr/bin/env python3
import os
from ase.io import read
from parameters import supercell_dim, sigma_phonon
from ...util.util import shell

def get_pdos():
    # PDOS
    config = read('POSCAR')
    n_atoms = len(config)  # n: number of
    print("chemical formula: "
          f"{config.get_chemical_formula(mode='reduce')}")
    # SiO2
    if config.get_chemical_formula(mode='hill', empirical=True) == 'O2Si':
        n_Si = int(n_atoms / 3)
        n_O = 2 * n_Si
        pdos_string = ''
        for i_atom in range(1, n_Si+1):  # 1 to n_Si
            pdos_string += f'{i_atom} '
        pdos_string += ', '
        for i_atom in range(n_Si+1, n_atoms+1):  # n_Si+1 to n_atoms
            pdos_string += f'{i_atom} '
    else:
        symbols = config.get_chemical_symbols()
        # print(f'{symbols =}')
        pdos_string = ''
        for index, symbol in enumerate(symbols):
            i_phonopy = index + 1  # start from 1
            if index >= 1 and symbols[index-1] != symbol:
                pdos_string += f', {i_phonopy} '
            else:
                pdos_string += f'{i_phonopy} '
    print(f'{pdos_string =}')

    shell(f'phonopy --mesh="1 1 1" --sigma="{sigma_phonon}" --hdf5 --readfc -p -s --dim="{supercell_dim}" --pdos="{pdos_string}" ')

    # Ref.
    # https://github.com/phonopy/phono3py/tree/master/example/Si-PBEsol


get_pdos()

shell(f'phonopy --mesh="1 1 1" --sigma="{sigma_phonon}" --hdf5 --readfc -p -s --dim="{supercell_dim}" ')

print('File written: partial_dos.pdf, projected_dos.dat')
print('  total_dos.pdf, total_dos.dat')


