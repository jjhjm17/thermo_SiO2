#!/usr/bin/env python3
"""This file reads distances to nearest atoms from
original_DFT_mediumCell_largest_alat_folder/POSCAR_supercell.vasp.
The 3rd (fcc, bcc) or 4th (hcp)  nearest neighbor is used to determine cutoff of MTP in
c_defineMore.py.

Restriction: periodic boundary condition is not considered.
"""

from ase.io import read
import numpy as np
# from ase import neighborlist
# from ase.neighborlist import NeighborList

# For quick MD
# mediumCell_largest_alat_POSCAR = \
#     './original_DFT_largest_alat_folder/POSCAR'
# For the medium supercell
mediumCell_largest_alat_POSCAR = \
    './original_DFT_largest_alat_folder/POSCAR_supercell.vasp'


# Ref:
# https://chemistry.stackexchange.com/questions/65395/finding-the-nearest-neighbour-from-an-atom-within-a-crystal
# get neighbors of atom 0

atoms = read(mediumCell_largest_alat_POSCAR)
distances = []
for i in range(len(atoms)):
    for j in range(i, len(atoms)):
        r1 = atoms.get_positions()[i]
        r2 = atoms.get_positions()[j]
        distance = np.sqrt( np.inner(r1 - r2, r1 - r2))
        # print('distance =', distance)
        distances.append(distance)

distances = np.unique( np.array(distances).round(decimals=4))
# distances = np.unique( np.array(distances).round(decimals=4) ).sort()
print('Neighbor distances from the nearest one: ', distances[:6], 'Ang')
