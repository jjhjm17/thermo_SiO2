#!python
"""
Functionalize a surface by inserting O–H groups on selected surface Si/Al atoms.
Breaks Si/Al–O–Si/Al bridges and replaces them with two OH groups pointing
inward toward the pore center (cylindrical mode) or slab mid-plane (slab mode).

Requirements:
- Surface atom indices in files `surface_Si_Al` and `surface_O`.
- Input structure in LAMMPS dump or data format.
- `specorder`: ['Si', 'O', 'Al'] for dump format.
- Output functionalized.data has atom order: ['Si', 'O', 'H', 'Al'].
"""

import argparse
import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.geometry import find_mic
from ase.geometry.geometry import naive_find_mic
from ase.neighborlist import NeighborList

def main():
    # -------------------------------
    # Argument parsing
    # -------------------------------
    parser = argparse.ArgumentParser(allow_abbrev=False, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("--dump", type=str, default="final_melt_quenched_structure.dump",
                        help="Input pore structure file")
    parser.add_argument("--format", type=str, default="lammps-dump-text",
                        help="ASE format for reading input")
    parser.add_argument("--seed", type=int, default=111,
                        help="Random seed for atom selection")
    parser.add_argument("--dist_O_H", type=float, nargs=2, default=[1.7, 1.0],
                        help="Distances (Å) for O insertion from Si/Al and H insertion from O")
    parser.add_argument("--safe_dist", type=float, default=1.5,
                        help="Minimum allowed distance (Å) from new atoms to existing atoms")
    parser.add_argument("--vac_shape", type=str, choices=["cylinder", "slab"], default="cylinder",
                        help="Vacuum shape for inward direction calculation")
    parser.add_argument("--vac_thick", type=float, default=15.0,
                        help="Vacuum thickness in Å (only for slab mode)")

    args = parser.parse_args()

    # -------------------------------
    # Store args in variables
    # -------------------------------
    safe_dist = args.safe_dist
    dist_O, dist_H = args.dist_O_H
    seed = args.seed
    dump_file = args.dump
    vac_shape = args.vac_shape
    vac_thick = args.vac_thick

    # -------------------------------
    # -----------------------------
    # --- helper function for MIC ---
    # -----------------------------
    def find_mic_custom(vec, cell, pbc):
        """
        If vac_shape = cylinder,
        Fast minimum image convention for orthorhombic cells. Otherwise
        use  ase's find_mic.
        vec: array-like, shape (3,) or (N,3)
        cell: ASE cell array
        Returns diffvec with same shape as vec
        It can be also probably used for hexagonal cells because the pore
        surface is along the z direction, and so only the z component
        matters at the surface.
        """
        if vac_shape == 'cylinder':
            return find_mic(vec, cell, pbc)
        else:  # slab, orthorhombic
            vec = np.asarray(vec)
            vec = np.atleast_2d(vec)
            diffvec, vlen = naive_find_mic(vec, cell=cell)
            diffvec = diffvec[0]  # We assume vec is a single vector, not vectors.
            return diffvec, vlen

    # -------------------------------
    # Read input structure
    # -------------------------------
    symbols = ['Si', 'O', 'Al']  # order in input dump
    sym_func = ['Si', 'O', 'H', 'Al']  # order for output

    if args.format == 'lammps-dump-text':
        atoms = read(dump_file, format='lammps-dump-text', specorder=symbols)
    elif args.format == 'lammps-data':
        Z_of_type = {1: 14, 2: 8, 3: 13}  # Si O Al
        atoms = read(dump_file, format='lammps-data', Z_of_type=Z_of_type)
        atoms.wrap()

    cell = atoms.get_cell()
    center_xy = (cell[0] + cell[1]) / 2  # for cylinder mode

    # -------------------------------
    # Additional slab parameters & checks
    # -------------------------------
    if vac_shape == "slab":
        # Rectangular cell check (orthogonal lattice)
        if not (np.allclose(cell[0, 1:], 0, atol=1e-8) and
                np.allclose(cell[1, [0, 2]], 0, atol=1e-8) and
                np.allclose(cell[2, :2], 0, atol=1e-8)):
            raise ValueError("Error: For vac_shape='slab', the simulation cell must be rectangular!")

        cell_z_len = cell[2][2]
        center_z_vac = cell_z_len - vac_thick / 2  # mid-plane of vacuum region
        slab_half_thickness = (cell_z_len - vac_thick) / 2



    # -------------------------------
    # Function to attempt OH insertion
    # -------------------------------

    # Precompute cutoff radius
    cutoff_radius = dist_O + dist_H + safe_dist + 0.01

    # Initialize NeighborList for the existing atoms
    cutoffs = [cutoff_radius / 2] * len(atoms)  # ASE uses half the cutoff for some internal reasons
    nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
    nl.update(atoms)

    def try_insert(i0, i1, shift=np.array(4 * [3 * [0.]])):
        """
        Attempt to insert OH groups between atoms i0 and i1.
        shift: offset vector (Å) to avoid overlaps.
        Returns:
          - If overlap detected: shift vector to adjust placement
          - If successful: (Atoms(O2), Atoms(H2))
        """
        p0 = atoms[i0].position
        p1 = atoms[i1].position

        if vac_shape == "cylinder":
            # Vector from atom to cylinder axis center (in XY plane only)
            p0_to_c = center_xy - p0
            p1_to_c = center_xy - p1
            p0_to_c[2] = 0
            p1_to_c[2] = 0
            p0_to_c /= np.linalg.norm(p0_to_c)
            p1_to_c /= np.linalg.norm(p1_to_c)

        else:  # slab mode with PBC-aware inward direction
            cell_z_len = cell[2][2]
            z_center = center_z_vac

            # For atom 0: minimal-image z displacement to mid-plane
            dz0 = z_center - p0[2]
            dz0 -= np.rint(dz0 / cell_z_len) * cell_z_len
            p0_to_c = np.array([0, 0, dz0])

            # For atom 1
            dz1 = z_center - p1[2]
            dz1 -= np.rint(dz1 / cell_z_len) * cell_z_len
            p1_to_c = np.array([0, 0, dz1])

            # Normalize → always points inward (shortest PBC path)
            p0_to_c /= np.linalg.norm(p0_to_c)
            p1_to_c /= np.linalg.norm(p1_to_c)

        # Place O and H positions
        pO0 = p0 + dist_O * p0_to_c + shift[0]
        pO1 = p1 + dist_O * p1_to_c + shift[1]
        pH0 = pO0 + dist_H * p0_to_c + shift[2]
        pH1 = pO1 + dist_H * p1_to_c + shift[3]

        # Check if any new atom is too close to existing atoms
        # for atom in atoms:
        #     pos = atom.position
        #     for i, p in enumerate([pO0, pO1, pH0, pH1]):
        #         diffvec, _ = find_mic_ortho(p - pos, cell=atoms.get_cell())
        #         diff = np.linalg.norm(diffvec)
        #         if diff < safe_dist:
        #             vec = np.array(4 * [3 * [0.]])
        #             v = 1.01 * (safe_dist - diff) * diffvec / diff
        #             vec[i] += v
        #             return vec

        # Check if any new atom is too close to existing atoms
        new_positions = np.array([pO0, pO1, pH0, pH1])
        for idx_new, pos_new in enumerate(new_positions):
            # Check neighbors of i0 and i1 first (most likely close)
            neighbors_i0, offsets_i0 = nl.get_neighbors(i0)
            neighbors_i1, offsets_i1 = nl.get_neighbors(i1)
            neighbors = np.unique(np.concatenate([neighbors_i0, neighbors_i1]))
            
            for neighbor in neighbors:
                neighbor_pos = atoms[neighbor].position
                diffvec, _ = find_mic_custom(pos_new - neighbor_pos, cell=atoms.get_cell(),
                                             pbc=atoms.get_pbc())
                difflen = np.linalg.norm(diffvec)
                if difflen < safe_dist:
                    vec = np.array(4 * [3 * [0.]])
                    v = 1.01 * (safe_dist - difflen) * diffvec / difflen
                    vec[idx_new] += v
                    return vec

        # No overlaps → return OH atoms
        aO2 = Atoms('O2', positions=[pO0, pO1])
        aH2 = Atoms('H2', positions=[pH0, pH1])
        return aO2, aH2

    # -------------------------------
    # Read surface atom lists
    # -------------------------------
    infile1 = 'surface_Si_Al'
    infile2 = 'surface_O'

    with open(infile1, 'r') as f:
        ind = [int(line.strip()) for line in f]

    with open(infile2, 'r') as f:
        O_neigh = [[int(x) for x in line.strip().split()] for line in f]

    # -------------------------------
    # Main functionalization loop
    # -------------------------------
    rng = np.random.default_rng(seed)
    checked = np.array([False for _ in range(len(ind))])
    to_del = []  # O atoms to remove
    to_add = []  # OH atoms to add
    ntries = 4


    while not all(checked):
        print(np.count_nonzero(checked == False))
        false_indices = np.where(checked == False)[0]
        rind = rng.choice(false_indices)  # pick random surface Si/Al
        checked[rind] = True
        Oind = O_neigh[rind]
        i0 = ind[rind]

        # Check all O neighbors for possible bridge breaking
        checked_O = np.array([False for _ in range(len(Oind))])
        while not all(checked_O):
            false_indices_O = np.where(checked_O == False)[0]
            rO = rng.choice(false_indices_O)
            checked_O[rO] = True
            iO = Oind[rO]

            # Look for another Si/Al sharing this O
            for j in range(len(ind)):
                if not checked[j] and iO in O_neigh[j]:
                    i1 = ind[j]
                    shift = np.array(4 * [3 * [0.]])
                    t = try_insert(i0, i1, shift=shift)
                    count = 0
                    while len(t) == 4 and count < ntries:
                        print("  ", count)
                        shift += t
                        t = try_insert(i0, i1, shift=shift)
                        count += 1
                    if len(t) == 2:
                        aO2, aH2 = t
                        # diffvec1, _ = find_mic(aO2[0].position - aH2[0].position,
                        #                        cell=atoms.get_cell(), pbc=atoms.get_pbc())
                        diffvec1, _ = find_mic_custom(aO2[0].position - aH2[0].position,
                                               cell=atoms.get_cell(), pbc=atoms.get_pbc())
                        # diffvec2, _ = find_mic(aO2[1].position - aH2[1].position,
                        #                        cell=atoms.get_cell(), pbc=atoms.get_pbc())
                        diffvec2, _ = find_mic_custom(aO2[1].position - aH2[1].position,
                                               cell=atoms.get_cell(), pbc=atoms.get_pbc())
                        diff1 = np.linalg.norm(diffvec1)
                        diff2 = np.linalg.norm(diffvec2)
                        if not diff1 < 0.99 * dist_H and not diff2 < 0.99 * dist_H:
                            checked_O = [True for _ in range(len(Oind))]
                            to_add.append(aO2)
                            to_add.append(aH2)
                            to_del.append(iO)
                            checked[j] = True
                            break

    # -------------------------------
    # Finalize structure and write outputs
    # -------------------------------
    to_keep = [i for i in range(len(atoms)) if i not in to_del]
    functionalized = atoms[to_keep]
    for a in to_add:
        functionalized += a

    # Sort atoms for output
    sort_indices = np.argsort([sym_func.index(sym) for sym in functionalized.get_chemical_symbols()])
    functionalized_sorted = functionalized[sort_indices]

    write('functionalized.xyz', functionalized_sorted, format='extxyz')
    write('functionalized.data', functionalized_sorted, format='lammps-data', specorder=sym_func)

    print("Functionalization complete.")


if __name__ == "__main__":
    main()
