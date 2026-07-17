#!/usr/bin/env python
"""This script makes lattice and parameters for the initial configration of Al-SiO2."""
import sys
import numpy as np
from shutil import copy
from ase import Atoms, Atom
from ase.io import write
from ase.build.attach import attach_randomly, attach
from a_parameters import cell_shape, alat, formula_num_Si, formula_num_Al, hole_radius as hole_radius_global
import a_parameters


def get_scaled_cell():
    """This function prints the number of atoms and returns the cell size."""

    mass_Si = 28.0855  # getAtomicMass.sh Si
    mass_O = 15.9994
    mass_Al = 26.9815386

    if cell_shape == 'cubic':
        clat = alat
        sin_gamma = 1  # gamma: angle between a and b vectors
    else:
        clat = a_parameters.clat
        sin_gamma = np.sqrt(3)/2

    print(
        f"input: {cell_shape = }, {alat = } Ang, {clat = } Ang, "
        f" {formula_num_Si = }, {formula_num_Al = }"
    )
    if formula_num_Al % 2 == 1:
        print(
            "Error: formula_num_Al should be a multiple of 2 because " "we add Al2O3."
        )
        sys.exit()
    formula_num_O = int(formula_num_Si * 2 + formula_num_Al / 2 * 3)
    print(f"\n output: {formula_num_O = }")

    mole_ratio_SiO2 = formula_num_Si / (formula_num_Si + formula_num_Al / 2)
    mole_ratio_Al2O3 = formula_num_Al / 2 / (formula_num_Si + formula_num_Al / 2)
    # Al_Si_ratio = formula_num_Al / formula_num_Si
    Al_Si_ratio = np.float64(formula_num_Al) / formula_num_Si
    # np.float64: gives inf when divided by 0
    print(f"{Al_Si_ratio = :.3f}  {mole_ratio_SiO2  = :.3f}")

    print(f"{mole_ratio_Al2O3 = :0.4f}")
    if formula_num_Si == 0: # Al2O3
        density_Al2O3_SiO2 = 3.97  # g/cm^3

        # Ref
        # https://pubchem.ncbi.nlm.nih.gov/compound/Aluminum-Oxide#section=Density
        # Soft white powder; transforms to corundum at 1200 °C; density: 3.97 g/cu cm; insoluble in water; soluble in acid; slightly soluble in alkaline solutions /Gamma-alumina/
        # Lide, D.R. CRC Handbook of Chemistry and Physics 88TH Edition 2007-2008. CRC Press, Taylor & Francis, Boca Raton, FL 2007, p. 4-45

        # Note the extrapolated value from the below formula gives a different
        # value, density_Al2O3_SiO2 = 3.120 g/cm^3

    else:  # contains SiO2
        density_SiO2 = 2.20  # g/cm^3
        density_with_Al2O3 = 2.43
        # density_with_Al2O3: exp. density of 25 mole % Al2O3, 2.43 g/cm^3
        # (SiO2)_(1-x) (Al2O3)_x, x=0.25
        # ref.  Structure of SiO2–Al2O3 glasses: Combined X-ray diffraction,
        # IR and Raman studies
        # Journal of Non-Crystalline Solids 351 (2005) 1032–1038
        # https://doi.org/10.1016/j.jnoncrysol.2005.01.014
        #
        # extrapolate
        # 2.20 + (2.43-2.20) * (100/25) = 3.12 g/cm^3
        # Al2O3 density: 3.99 g/cm^3 (Wikipedia)
        # The difference probably comes from the different coordination numbers.
        density_Al2O3_SiO2 = (
            density_SiO2 + (density_with_Al2O3 - density_SiO2) * mole_ratio_Al2O3 / 0.25
        )

    # volume = alat**2 * clat * np.sqrt(3) / 2  # Ang^3
    hole_radius = hole_radius_global  # make local for locals() later
    volume = alat**2 * clat * sin_gamma - np.pi * hole_radius**2 * clat   # Ang^3  out of  cylinder
    print(f"{density_Al2O3_SiO2 = :0.3f} g/cm^3   {volume = :.3f} Ang^3")

    # gnu units
    # You have: amu/angstrom^3
    # You want: g/cm^3
    #         * 1.6605391
    amu_per_angstrom_3__to__g_per_cm_3 = 1.6605391
    density_amu_per_angstrom_3 = density_Al2O3_SiO2 / amu_per_angstrom_3__to__g_per_cm_3
    mass_cell = density_amu_per_angstrom_3 * volume

    # mass_cell should be a multiple of the mass of chemical formula
    mass_formula = (
        mass_Si * formula_num_Si + mass_O * formula_num_O + mass_Al * formula_num_Al
    )
    num_formula_float = mass_cell / mass_formula
    num_formula = int(np.round(mass_cell / mass_formula))
    # number of atoms of the chemical formula in the cell.
    num_Si = num_formula * formula_num_Si
    num_SiO2 = num_Si
    num_O = num_formula * formula_num_O
    num_Al = num_formula * formula_num_Al
    num_Al2O3 = int(num_Al/2)
    num_atom = num_Si + num_O + num_Al
    num_atom_dict = dict(Si=num_Si, O=num_O, Al=num_Al)

    print(f'{hole_radius = }, {volume = }')
    print(f"{mass_formula = :.2f}, {mass_cell = :.2f}, {num_formula = }")
    print(f"{num_Si = }, {num_O = }, {num_Al = }, {num_atom = }")
    print(f"{num_SiO2 = }, {num_Al2O3 = }")

    # Adjust the cell size to obtain the experimental density
    # rho = m / V
    # volume * scale**3 = volume_scaled
    # num_formula_float / volume = num_formula / volume_scaled
    # scale = (volume_scaled / volume)**(1/3) = (num_formula / num_formula_float)**(1/3)
    scale = (num_formula / num_formula_float)**(1/3)
    print("\n We adjust the cell size by 'scale' to obtain the experimental density.")
    alat_scaled = alat * scale
    clat_scaled = clat * scale
    print(f'{scale = :.4f}')
    print(f'{alat_scaled = :.4f} Ang, {clat_scaled = :.4f} Ang.')

    # # mean distance for cutoff of create_random
    # # V = 4/3 pi r^2    r = (3/(4 pi) V )^(1/3)
    # mean_dist_O = (3 / (4 * np.pi) * volume / num_O)**(1/3)
    # print(f'{mean_dist_O = :.4f} Ang')


    return alat_scaled, clat_scaled, num_atom_dict, scale, cell_shape


def write_lammps_lattice_dataf(alat, clat, cell_shape):
    """This function writes the lammps datafile of a hexagonal or cubic cell."""
    if cell_shape == 'hexagonal':
        with open('lattice.dataf', 'w') as fout:
            fout.write('Hexagonal cell lattice\n')
            fout.write('\n')
            fout.write('3 atom types\n')
            fout.write(f'0 {alat}  xlo xhi\n')
            fout.write(f'0 {np.sqrt(3)/2 * alat}  ylo yhi\n')
            fout.write(f'0 {clat}  zlo zhi\n')
            fout.write(f'{alat/2} 0 0  xy xz yz\n')
            fout.write('\n')
            fout.write('\n')
            fout.write('Atoms\n')
            # fout.write('\n')
            # fout.write('\n')
            # fout.write('\n')
            # fout.write('Atoms\n')
    elif cell_shape == 'cubic':
        with open('lattice.dataf', 'w') as fout:
            fout.write('Cubic cell lattice\n')
            fout.write('\n')
            fout.write('3 atom types\n')
            fout.write(f'0 {alat}  xlo xhi\n')
            fout.write(f'0 {alat}  ylo yhi\n')
            fout.write(f'0 {alat}  zlo zhi\n')
            fout.write(f'0 0 0  xy xz yz\n')
            fout.write('\n')
            fout.write('\n')
            fout.write('Atoms\n')
            # fout.write('\n')
            # fout.write('\n')
            # fout.write('\n')
            # fout.write('Atoms\n')
    else:
        print('Error: unknown cell_shape.')
        sys.exit()

    print('written: lattice.dataf')



# def make_random_initial(alat, clat, num_atom_dict):
#     """This script makes random initials with cutoff."""
# def make_lattice(alat, clat):
#     """This script makes random initials with cutoff."""
#     # cutoff = 2  # Ang
#     # print(f'\n make_random_initial')
#     config = Atoms(pbc=True)
#     config.cell = [alat, alat, clat, 90, 90, 60]
#     print(f'{config.get_cell() = }')
#     # atoms = ['Si', 'O', 'Al']
#     # # class rng_fix(np.random.default_rng):
#     # #     def rand(self):
#     # #         self.random()
#     # rng = np.random
#     # rng.default_rng(123456)
#     # # rng = rng_fix(123456)
#     # for atom in atoms:
#     #     for ind in range(num_atom_dict[atom]):
#     #         if atom == atoms[0] and ind == 0:
#     #             config.append(Atom(atom))
#     #         else:
#     #             config = attach_randomly(config, Atoms(atom, cell=config.cell,
#     #                                                    pbc=True), cutoff, rng)
# 
#     write('lattice.dataf', config, format='lammps-data')
#     print('\n written: lattice.dataf')
#     print('Manually change the number of atom types from 0 to 3 for the next step. and change xy to the positive value')

def round_input(num):
    return np.round(num, 4)


def write_lammps_input(alat, clat, num_atom_dict, scale):
    """This function reads input_file_blank file (default in.file.blank) and writes in.file.
    For example, for input_file_blank = 'in.header.blank',
    it writes in.header file."""
    variables = ['num_SiO2', 'num_Al2O3', 'hole_center_x', 'hole_center_y', 'hole_radius', 'alat_start', 'clat_start', 'xy_start']
    num_SiO2 = num_atom_dict['Si']
    # num_O = num_atom_dict['O']
    # num_Al = num_atom_dict['Al']
    num_Al2O3 = int(num_atom_dict['Al'] / 2)

    alat_start = round_input(alat)
    clat_start = round_input(clat)
    if cell_shape == 'hexagonal':
        hole_center_x = (alat + alat/2  ) / 2  # (vec_a + vec_b)/2
        hole_center_y = alat * np.sqrt(3) / 2 / 2
        xy_start = alat / 2
    elif cell_shape == 'cubic':
        hole_center_x = alat / 2
        hole_center_y = alat / 2
        xy_start = 0
    else:
        print('Error: unknown cell_shape.')
        sys.exit()
    hole_radius = hole_radius_global * scale  # make a local variable for locals()

    hole_center_x = round_input(hole_center_x)
    hole_center_y = round_input(hole_center_y)
    xy_start = round_input(xy_start)
    hole_radius = round_input(hole_radius)

    if hole_radius < 1e-6:
        hole_radius = 1e-6  # Otherwise lammps cylinder command makes error.

    if hasattr(a_parameters, 'input_file_blank'):
        input_file_blank = a_parameters.input_file_blank
        input_file = input_file_blank.replace('.blank', '')
    else:
        input_file_blank = 'in.file.blank'
        input_file = 'in.file'
    if hasattr(a_parameters, 'out_folder'):
        out_folder = a_parameters.out_folder
    else:
        out_folder = '.'
    with open(f'./{input_file_blank}', 'r') as fin:
        with open(f'{out_folder}/{input_file}', 'w') as fout:
            for line in fin.readlines():
                for variable in variables:
                    if f'xxx_{variable}_xxx' in line:
                        line = line.replace(f'xxx_{variable}_xxx', f'{locals()[variable]}')
                fout.write(line)
    print(f'\n read: {input_file_blank}')
    print(f' written: {out_folder}/{input_file}')


def make_initial_config():
   alat_scaled, clat_scaled, num_atom_dict, scale, cell_shape = get_scaled_cell()
   # make_random_initial(alat_scaled, clat_scaled, num_atom_dict)
   # write_lammps_lattice_dataf(alat_scaled, clat_scaled, cell_shape)
   write_lammps_input(alat_scaled, clat_scaled, num_atom_dict, scale)


if __name__ == "__main__":
    make_initial_config()
