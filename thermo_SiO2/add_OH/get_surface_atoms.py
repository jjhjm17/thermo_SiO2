#!python
import argparse
import itertools
import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.geometry.geometry import complete_cell, minkowski_reduce, wrap_positions, Cell, pbc2pbc
from ase.neighborlist import NeighborList
# from ase.geometry import find_mic


def general_find_mic_fast(v, cell, rcell, pbc=True):
    """Finds the minimum-image representation of vector(s) v. Using the
    Minkowski reduction the algorithm is relatively slow but safe for any cell.
    rcell: reduced cell
    """

    cell = complete_cell(cell)
    # rcell, _ = minkowski_reduce(cell, pbc=pbc)
    positions = wrap_positions(v, rcell, pbc=pbc, eps=0)

    # In a Minkowski-reduced cell we only need to test nearest neighbors,
    # or "Voronoi-relevant" vectors. These are a subset of combinations of
    # [-1, 0, 1] of the reduced cell vectors.

    # Define ranges [-1, 0, 1] for periodic directions and [0] for aperiodic
    # directions.
    ranges = [np.arange(-1 * p, p + 1) for p in pbc]

    # Get Voronoi-relevant vectors.
    # Pre-pend (0, 0, 0) to resolve issue #772
    hkls = np.array([(0, 0, 0)] + list(itertools.product(*ranges)))
    vrvecs = hkls @ rcell

    # Map positions into neighbouring cells.
    x = positions + vrvecs[:, None]

    # Find minimum images
    lengths = np.linalg.norm(x, axis=2)
    indices = np.argmin(lengths, axis=0)
    vmin = x[indices, np.arange(len(positions)), :]
    vlen = lengths[indices, np.arange(len(positions))]
    return vmin, vlen

def find_mic_fast(v, cell, rcell, pbc=True):
    """Finds the minimum-image representation of vector(s) v using either one
    of two find mic algorithms depending on the given cell, v and pbc."""

    cell = Cell(cell)
    pbc = cell.any(1) & pbc2pbc(pbc)
    dim = np.sum(pbc)
    v = np.asarray(v)
    single = v.ndim == 1
    v = np.atleast_2d(v)

    vmin, vlen = general_find_mic_fast(v, cell, rcell, pbc=pbc)

    if single:
        return vmin[0], vlen[0]
    else:
        return vmin, vlen

def main():
    # parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(allow_abbrev=False, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        "--dump",
        type = str,
        default = "final_melt_quenched_structure.dump",
        help = "pore structure, lammps-dump-text format, specorder: ['Si', 'O', 'Al']"
    )

    parser.add_argument(
        "--format",
        type = str,
        default = "lammps-dump-text",
        help = "ase.io.read format, specorder: ['Si', 'O', 'Al']"
    )

    parser.add_argument(
        "--dr",
        type = float,
        default = 3,
        help = "thickness of cylindrical sections (Angstrom) to include surface atoms"
    )

    parser.add_argument(
        "--dR",
        type = float,
        default = 3,
        help = "lateral extension of cylindrical sections (Angstrom) to include surface atoms"
    )

    parser.add_argument(
        "--cutoff",
        type = float,
        default = 2.1,
        help = "cutoff (Angstrom) for Si-O and Al-O bonds"
    )

    parser.add_argument(
        "--check",
        type = bool,
        default = False,
        help = "write surface.xyz with surface atoms for check in ovito"
    )

    parser.add_argument(
        "--symbols",
        type = str,
        default = 'Si O H Al',
        help = "Symbols, default: 'Si O H Al'"
    )

    args = parser.parse_args()
    dr = args.dr
    dR = args.dR
    cutoff = args.cutoff
    check_with_ovito = args.check
    dump_file = args.dump


    # symbols = ['Si', 'O', 'Al']
    symbols = args.symbols.split()
    outfile1 = 'surface_Si_Al'
    outfile2 = 'surface_O'

    # atoms = read(dump_file, format='lammps-dump-text', specorder=symbols)
    format = args.format  # Note that in pdb args.format does not work.
    if format == 'lammps-dump-text':
        atoms = read(dump_file, format=format, specorder=symbols)
    elif format == 'lammps-data':
        # Z_of_type = {1:14, 2:8, 3:13}  # Si O Al
        if symbols == ['Si', 'O', 'H', 'Al']:
            Z_of_type = {1:14, 2:8, 3:1, 4:13}  # Si O H Al
        atoms = read(dump_file, format=format, Z_of_type=Z_of_type)
        atoms.wrap()
    cell = atoms.get_cell()
    rcell, _ = minkowski_reduce(cell, pbc=atoms.get_pbc())
    center = (cell[0] + cell[1]) / 2

    ind = []
    dist = []

    for atom in atoms:
      if atom.symbol == 'Si' or atom.symbol == 'Al':
        i = atom.index
        p = atom.position
        r = p - center
        d = np.linalg.norm(r[:2])
        ind.append(i)
        dist.append(d)

    dist, ind = zip(*sorted(zip(dist,ind)))
    checked = [False for i in range(len(dist))]

    while not all(checked):
      index = next((i for i, val in enumerate(checked) if not val), None)
      print(index, len(checked))
      i0 = ind[index]
      d0 = dist[index]
      checked[index] = True
      p0 = atoms[i0].position
      r0 = p0 - center
      dalpha = dR / d0

      to_delete = []
      for i in range(1, len(dist)):
        i1 = ind[i]
        d1 = dist[i]
        p1 = atoms[i1].position
        r1 = p1 - center
        v = p1 - p0
        v_mic, _ = find_mic_fast(v, cell=atoms.get_cell(), rcell=rcell, pbc=atoms.get_pbc())
        dc = v_mic[2]
        
        r0_ = r0[:2]
        r1_ = r1[:2]
        n0_ = np.linalg.norm(r0_)
        n1_ = np.linalg.norm(r1_)

        # to avoid warning, we cut at 1; in rare cases argument can be numerically above 1
        d = np.min([np.dot(r0_,r1_) / (n0_*n1_), 1])
        alpha = np.arccos(d)
        
        if dc < dR and np.abs(alpha) < dalpha and np.abs(d1-d0) > dr:
          to_delete.append(i)
      
      ind = [x for i, x in enumerate(ind) if i not in to_delete]
      dist = [x for i, x in enumerate(dist) if i not in to_delete]
      checked = [x for i, x in enumerate(checked) if i not in to_delete]

    # 2. Finding oxygen neighbors of these surface Si/Al atoms
    #     O_neigh = [[] for i in range(len(ind))]
    #     for i in range(len(ind)):
    #       print(i)
    #       i0 = ind[i]
    #       for atom in atoms:
    #         if atom.symbol == "O":
    #           i1 = atom.index
    #           d = atoms.get_distance(i0, i1, mic=True)
    #           d = d_s[i1]
    #           if d < cutoff:
    #             O_neigh[i].append(i1)

    # 2. Finding oxygen neighbors of these surface Si/Al atoms
    # faster version
    # ref: chat gpt
    cutoff_half_bond = cutoff / 2
    cutoff_dict = {'Si': cutoff_half_bond, 'Al': cutoff_half_bond,
                    'O': cutoff_half_bond}  # per-type cutoffs if needed
    cutoffs = [cutoff_dict[sym] for sym in atoms.get_chemical_symbols()]

    # Build neighbor list (True = self-interaction excluded)
    nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
    nl.update(atoms)

    O_neigh = [[] for _ in range(len(ind))]
    for idx, si_al_index in enumerate(ind):
        indices, offsets = nl.get_neighbors(si_al_index)
        for j, offset in zip(indices, offsets):
            if atoms[j].symbol == "O":
                O_neigh[idx].append(j)


    if check_with_ovito:
      surind = np.array([x for l in O_neigh for x in l])
      surind = np.concatenate((np.array(ind), surind))
      surface = atoms[surind]
      write('surface.xyz', surface, format='extxyz')

    with open(outfile1, 'w') as f:
      for i in ind:
        f.write(str(i) + '\n')

    with open(outfile2, 'w') as f:
      for i in O_neigh:
        for j in i:
          f.write(str(j) + ' ')
        f.write('\n')


if __name__ == "__main__":
    main()
