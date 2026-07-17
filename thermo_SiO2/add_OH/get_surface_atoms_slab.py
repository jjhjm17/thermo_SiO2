#!python
"""
We assume that the cell is rectangular.

to_delete contains Si/Al atoms that are in nearly the same lateral direction and
height as the reference atom but are further or nearer in radial distance by
more than dr.  The goal is to avoid redundant atoms in the same vertical/angular
neighborhood, possibly to get only the outer surface Si/Al atoms.


to_delete is not the surface atoms — it’s the “discard list” of redundant Si/Al
atoms removed from consideration.

The atoms that survive after removing all to_delete entries are the core surface
Si/Al atoms, which are then combined with their oxygen neighbors to form
surface.
"""


import argparse
import numpy as np
from ase import Atoms
from ase.io import read, write
# from ase.geometry import find_mic
from ase.neighborlist import NeighborList

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

# parser.add_argument(
#     "--dr",
#     type = float,
#     default = 3,
#     help = "thickness of cylindrical sections (Angstrom) to include surface atoms"
# )

parser.add_argument(
    "--dz",
    type = float,
    default = 3,
    help = "thickness of vacuum regions (Angstrom) to include surface atoms"
)

parser.add_argument(
    "--vac_thick",
    type = float,
    default = 15,
    help = "thickness of vacuum region (Angstrom)"
)

parser.add_argument(
    "--dR",
    type = float,
    default = 3,
    # help = "lateral extension of cylindrical sections (Angstrom) to include surface atoms"
    help = "lateral extension of vacuum regions along x, y directions (Angstrom) to include surface atoms"
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
    "--move_vacuum_to_center",
    type = bool,
    default = False,
    help = "move the vacuum region to the center in the surface.xyz file."
)

args = parser.parse_args()
vac_thick = args.vac_thick
dz = args.dz
dR = args.dR
cutoff = args.cutoff
check_with_ovito = args.check
dump_file = args.dump


symbols = ['Si', 'O', 'Al']
outfile1 = 'surface_Si_Al'
outfile2 = 'surface_O'

# atoms = read(dump_file, format='lammps-dump-text', specorder=symbols)
format = args.format  # Note that in pdb args.format does not work.
if format == 'lammps-dump-text':
    atoms = read(dump_file, format=format, specorder=symbols)
elif format == 'lammps-data':
    Z_of_type = {1:14, 2:8, 3:13}  # Si O Al
    atoms = read(dump_file, format=format, Z_of_type=Z_of_type)
    atoms.wrap()
cell = atoms.get_cell()
# center = (cell[0] + cell[1]) / 2

# check if the simulation cell is rectangular by inspecting the off-diagonal terms of the cell matrix.
# ref: chat gpt
if not np.allclose(cell[0,1:], 0, atol=1e-8) or \
   not np.allclose(cell[1,[0,2]], 0, atol=1e-8) or \
   not np.allclose(cell[2,:2], 0, atol=1e-8):
    raise ValueError("Error: Simulation cell is not rectangular!")

# Define your vacuum center
center_z_vac = cell[2][2] - vac_thick / 2
slab_half_thickness = (cell[2][2] - vac_thick) / 2

# Build initial list of Si/Al atoms
ind = []
zdist = []
for atom in atoms:
  if atom.symbol == 'Si' or atom.symbol == 'Al':
    i = atom.index
    p = atom.position
    # r = p - center
    # d = np.linalg.norm(r[:2])
    # dz = abs(p[2] - center_z_vac)
    # d: distance from the center plane of the vacuum region
    ind.append(i)
    # dist.append(d)
    zdist.append(abs(p[2] - center_z_vac))

# Sort by vertical distance from the vacuum center
zdist, ind = zip(*sorted(zip(zdist,ind)))
checked = [False for _ in range(len(zdist))]

while not all(checked):
  # index = next((i for i, val in enumerate(checked) if not val), None)

  # Find the index of the first atom in checked that hasn't been marked as processed yet.
  # index: an atom in the list that is not checked
  index = next(i for i, val in enumerate(checked) if not val)
  print(index, len(checked))
  i0 = ind[index]
  z0 = atoms[i0].position[2]
  # d0 = dist[index]
  checked[index] = True
  p0 = atoms[i0].position
  # # r0 = p0 - center
  # # dz0 = p0[2] - center_z_vac
  # dalpha = dR / d0
  to_delete = []

  # for i in range(1, len(zdist)):
    # i1 = ind[i]
    # # d1 = dist[i]
    # p1 = atoms[i1].position
    # r1 = p1 - center
    # # d1 = p1[2] - center_z_vac
    # v = p1 - p0
    # v_mic, _ = find_mic(v, cell=atoms.get_cell(), pbc=atoms.get_pbc())
    # dc = v_mic[2]
    
    # r0_ = r0[:2]
    # # d0_ = d0[:2]
    # r1_ = r1[:2]
    # # d1_ = d1[:2]
    # n0_ = np.linalg.norm(r0_)
    # # n0_ = np.linalg.norm(d0_)
    # n1_ = np.linalg.norm(r1_)
    # # n1_ = np.linalg.norm(d1_)

    # # to avoid warning, we cut at 1; in rare cases argument can be numerically above 1
    # d = np.min([np.dot(r0_,r1_) / (n0_*n1_), 1])
    # alpha = np.arccos(d)
    
    # if dc < dR and np.abs(alpha) < dalpha and np.abs(d1-d0) > dr:
    #   to_delete.append(i)


  def find_mic_rect(v, cell, pbc):
    """Find the minimum image convention for a vector v with respect to a rectangular cell.
    This is faster than find_mic for rectangular cells."""
    if not pbc.any():
        return v, 0.0
    v_mic = np.zeros_like(v)
    for i in range(3):
        if pbc[i]:
            v_mic[i] = v[i] - np.round(v[i] / cell[i, i]) * cell[i, i]
        else:
            v_mic[i] = v[i]
    return v_mic, 0.0

  # find_mic_rect test
  def test_find_mic_rect():
    v=[0, 0, 12]
    cell=np.array([[10, 0, 0], [0, 10, 0], [0, 0, 10]])
    pbc=np.array([True, True, True])
    v_min_test = find_mic_rect(v, cell, pbc)
    print(f'v = {v}, cell = {cell}, pbc = {pbc} ')
    print('find_mic test  ,', v_min_test)

    v=[0, 0, -2]
    v_min_test = find_mic_rect(v, cell, pbc)
    print(f'v = {v}, cell = {cell}, pbc = {pbc} ')
    print('find_mic test  ,', v_min_test)

  # test_find_mic_rect()


  # Loop through all atoms and check their distance to the current atom
  # with index 'index'.
  for i in range(len(zdist)):
    if i == index:
      continue
    i1 = ind[i]  # index of the atom to check
    z1 = atoms[i1].position[2]  # z-coordinate of the atom to check

    p1 = atoms[i1].position
    v = p1 - p0  # vector
    # v_mic, _ = find_mic(v, cell=atoms.get_cell(), pbc=atoms.get_pbc())
    v_mic, _ = find_mic_rect(v, cell=atoms.get_cell(), pbc=atoms.get_pbc())
    # dx = abs(atoms[i1].position[0] - atoms[i0].position[0])
    # dy = abs(atoms[i1].position[1] - atoms[i0].position[1])

    # # apply periodic boundary conditions for x,y if needed:
    # dx = min(dx, cell[0,0] - dx)
    # dy = min(dy, cell[1,1] - dy)
    d_xy = np.linalg.norm(v_mic[:2])
    dist_z_vec = abs(v_mic[2])

    # deletion rule for rectangular region:
    # - vertical distance small (same horizontal "layer")
    # - lateral distance small
    # - big difference in zdist from vacuum plane
    # if dx < dR and dy < dR and abs(zdist[i] - zdist[index]) > dz:
    if d_xy < dR and dz < dist_z_vec < slab_half_thickness:
        to_delete.append(i)    
  
  ind = [x for i, x in enumerate(ind) if i not in to_delete]
  zdist = [x for i, x in enumerate(zdist) if i not in to_delete]
  checked = [x for i, x in enumerate(checked) if i not in to_delete]


# 2. Finding oxygen neighbors of these surface Si/Al atoms
# O_neigh = [[] for i in range(len(ind))]
# for i in range(len(ind)):
#   print(i)
#   i0 = ind[i]
#   for atom in atoms:
#     if atom.symbol == "O":
#       i1 = atom.index
#       d = atoms.get_distance(i0, i1, mic=True)
#       if d < cutoff:
#         O_neigh[i].append(i1)

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
  if args.move_vacuum_to_center:
    # Move the vacuum region to the center in the surface.xyz file.
    # This is useful for visualization in ovito.
    surface.translate([0, 0, - (center_z_vac - cell[2, 2] / 2)])
    surface.wrap()
  write('surface.xyz', surface, format='extxyz')

with open(outfile1, 'w') as f:
  for i in ind:
    f.write(str(i) + '\n')

with open(outfile2, 'w') as f:
  for i in O_neigh:
    for j in i:
      f.write(str(j) + ' ')
    f.write('\n')


