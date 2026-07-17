# alat, clat, hole_radius given here are approximate. Will be rescaled
# due to finite number of atoms to match the experimental density.

cell_shape = 'cubic'
# cell_shape = 'hexagonal'

# alat = 35.4  # Angst, 3000 atoms
# alat = 22  # Angst, 700 atoms
# Ref borosilicate: 3000 atoms, NVT exp density (2 commercial glasses).
# Lee, Kuo-Hao, ..., John C. Mauro et al. "Evaluation of classical interatomic potentials for molecular dynamics simulations of borosilicate glasses." Journal of Non-Crystalline Solids 528 (2020): 119736. 
# https://doi.org/10.1016/j.jnoncrysol.2019.119736

#
# alat = 30  # Angst
# clat = 15  # Angst
# alat = 11  # Angst  small for active learning, vasp, larger than cutoff 5 Ang
# alat = 12  # Angst  small for active learning, vasp, larger than cutoff 5 Ang
# alat = 14  # Angst  small for active learning, vasp, aim for 200 atoms
alat = 14.5  # Angst  small for active learning, vasp, aim for 200 atoms
# clat = 15  # Angst

# clat = 11  # Angst  small for active learning, vasp, larger than cutoff 5 Ang

# hole_radius = 9  # Angst
# 30 / 2 / 2  + 1.5 = 9
# when we functionalize with OH, the hole will become smaller, so add 2.
# the size of the hole is similar to the SiO2 wall thickness.
# Si-O bond length: 1.6 (Google), O-H: 1, Si-O-H: 2.6, increase to account
# the radius

# hole_radius = 4  # Angst
hole_radius = 1e-6  # Angst, 0 gives error.
# 15 (alat) / 2 / 2  + 1.5 =~ 5, cell is small, 15 (alat) / 2 / 2 = 3.75 ~ 4


# Al_Si_ratio = 1/6.  # 1 : 6 as in the experiment for non-porous material
# formula_num_Si = 12  # number of Si in the reduced chemical formula
# formula_num_Si = 100 * 2  # number of Si in the reduced chemical formula
formula_num_Si = 2  # number of Si in the reduced chemical formula
# reduced: Si4O2 -> SiO2
# formula_num_Al = 2  # multiple of 2
# formula_num_Al = 17 * 2  # multiple of 2, consistent with Tom's 700 atoms
formula_num_Al = 0  # 
# formula_num_O = formula_num_Si * 2 + formula_num_Al / 2 * 3
# 10 * 2 + 2 / 2 * 3 = 20 + 3 = 23

input_file_blank = 'in.header.blank'
out_folder = '../c_melt'
