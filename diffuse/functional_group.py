"""This script counts the OH functional groups."""

from typing import List
from collections import Counter  # multiset
import numpy as np
from ase import Atoms
from ase.geometry.analysis import Analysis
from ase.neighborlist import NeighborList
from ase.neighborlist import neighbor_list
# from ase.neighborlist import build_neighbor_list
from thermo_SiO2.mlip.read import read_SiO2


def make_cutoffs(config: Atoms) -> List[float]:
    """This function makes the cutoff for the neighbor list."""
    shannon_radius = {'H': -0.38, 'O': 1.35, 'Si': 0.26, 'Al': 0.48}
        # http://abulafia.mt.ic.ac.uk/shannon/ptable.php
        # coordination number: H: 1, O: 2, Si: 4, Al: 4
    O2_bond_length = 1.2  # Ang
    reduce_O = shannon_radius['O'] - O2_bond_length/2
    # reduce_O = 0.4  # reduce O's radius so that H's radius becomes positive.
    positive_radius = {symbol: radius + reduce_O for (symbol, radius)
                           in shannon_radius.items() if symbol != 'O'}
    positive_radius['O'] = shannon_radius['O'] - reduce_O
    add = 0.15  # 2*add = 0.3 is larger than the half width of peak in RDF
    # made by ovito.
    positive_radius = {symbol: radius + add for symbol, radius in positive_radius.items()}
    cutoffs = [positive_radius[atom.symbol] for atom in config]

    # vdw_radius = {'H': 1.2, 'O': 1.52, 'Si': 2.1, 'Al': 2}
    # # ref:  ovito default
    # cutoffs = [vdw_radius[atom.symbol] for atom in config]
    return cutoffs


def make_all_bonds(config: Atoms) -> List[List[List[int]]]:
    """This function makes all_bonds inn config."""
    # i: first atom index, j: second atom index
    # nl_ij: [i's, j's]
    #       = [[bond_1_atom_index_1, bond_2_atom_index_1, ...],
    #         [bond_1_atom_index_2, bond_2_atom_index_2, ...]]
    # i's: ascending order

    # Make all_bonds from nl_ij,
    # because pairwise cutoff can be used only by neighbor_list funciton,
    # which makes nl_ij only and not all_bonds.
    #
    # all_bonds = [[[atom_1_bonded_atom_1, atom_1_bonded_atom_2, ...],
    #               [atom_2_bonded_atom_1, atom_2_bonded_atom_2, ...],
    #               ...], ]
    # all_bonds[0]: 1st snapshot

    nl_ij = neighbor_list('ij', config, cutoff=2.5)

    all_bonds = [[] for _ in range(len(config))]
    all_bonds = [all_bonds]  # [1st snapshot]
    for i_nl in range(len(nl_ij[0])):
        i_atom = nl_ij[0][i_nl]
        j_atom = nl_ij[1][i_nl]
        all_bonds[0][i_atom].append(int(j_atom))
        all_bonds[0][j_atom].append(int(i_atom))
    return all_bonds


def count_groups(file: str, atom_symbols: str, index: str = ":"):
    """This function counts Si-OH-Al groups (for now)."""
    configs = read_SiO2(file, atom_symbols, index)

    n_Al_OH_Si_s = []
    n_H3O_s = []
    VERBOSE = True
    # VERBOSE = False
    cutoffs = make_cutoffs(configs[0])
    for config in configs:
        O_indexes = []
        H3O_indexes = []

        nl = NeighborList(cutoffs, skin=0.1)
        nl.update(config)
        # neighbor_list = build_neighbor_list(config, cutoffs=cutoffs)
        ana = Analysis(configs, nl=nl)
        all_bonds = ana.all_bonds
        symbols = np.array(config.get_chemical_symbols())
        pos = config.get_positions()

        # all_bonds = make_all_bonds(config)


        for i_atom in range(len(config)):
            if symbols[i_atom] == "O":
                bonded_atom_i_s = all_bonds[0][i_atom]
                n_bonds = len(bonded_atom_i_s)
                if n_bonds >= 3:  # 3-bonded O
                    bonded_atoms = Counter(symbols[bonded_atom_i_s])
                    # if VERBOSE:
                    #     print(f'\n {i_atom = }')
                    #     print(all_bonds[0][i_atom])
                    print(bonded_atoms)
                    if n_bonds == 3:
                        if bonded_atoms == Counter(["H", "Al", "Si"]):
                            if VERBOSE:
                                print(f'{pos[i_atom] = }')
                            # print(pos[bonded_atom_i_s])
                            # print(bonded_atoms)
                            O_indexes.append(i_atom)
                        elif bonded_atoms == Counter(["H", "H", "H"]):
                            H3O_indexes.append(i_atom)
                    elif n_bonds > 3:
                        print('n_bonds > 3')
                        print(f'{n_bonds = }')
                        print(f'bonded_atoms = {bonded_atoms}')
        n_Al_OH_Si_s.append(len(O_indexes))
        n_H3O_s.append(len(H3O_indexes))
    # breakpoint()
    np.savetxt(
        "n_Al_OH_Si_s.txt",
        np.column_stack((np.linspace(0, len(configs) - 1, len(configs)), n_Al_OH_Si_s, n_H3O_s)),
        header="index   number of groups Al_OH_Si    H3O",
        fmt="%d %d %d",
    )

    return configs

