#!/usr/bin/env python3
"""This script makes folders for EOS fitting, for many volumes, to run vasp."""

import os
import sys
import shutil
import subprocess
import numpy as np
from ase.io import read, write
from .util import get_sample_folder_name
from ...util.util import check_link_symlink
from a_parameters import (calculation_folder as calc_folder,
                          unitcell_structure, use_vdw_kernel_file, num_samples,
                          min_vol_percent, max_vol_percent, num_samples,
                          POSCAR_prefix)

def get_alat_from_V(volume, unitcell_structure):
    """This function returns the lattice constant of a supercell from the
    volume."""
    if unitcell_structure == 'amorphous_cubic':
        alat = volume**(1/3)
    else:
        print("Error: for now only unitcell_structure = 'amorphous_cubic' is"
              " supported.")
        sys.exit()
    return alat


def get_V_from_alat(alat, unitcell_structure):
    """This function returns the volume of a supercell from the lattice
    constant."""
    if unitcell_structure == 'amorphous_cubic':
        vol = alat**3
    else:
        print("Error: for now only unitcell_structure = 'amorphous_cubic' is"
              " supported.")
        sys.exit()
    return vol


def get_alat_range(template_folder):
    """This function returns alat range for EOS."""
    # approx_eq_config = read('../template/POSCAR')
    # approx_eq_config = read(f'../template/{POSCAR_prefix}_0')
    approx_eq_config = read(f'{template_folder}/{POSCAR_prefix}_0')
    # cell volume
    vol_0 = approx_eq_config.get_volume() # approximate equilibrium volume
    # min_vol_percent = -8  # percent of volume
    # max_vol_percent = 12  # percent of volume

    # # Approximately consistent with Axel's one, SiO2 is very soft.
    # min_vol_percent = -12  # percent of volume
    # max_vol_percent = 18  # percent of volume
    vol_min_first  = vol_0 * (1 + min_vol_percent/100)
    vol_max_first  = vol_0 * (1 + max_vol_percent/100)
    print(f'volume {min_vol_percent = } %')
    print(f'volume {max_vol_percent = } %')
    print(f'{vol_0 = :.6g} Ang^3')
    print(f'{vol_min_first = :.6g} Ang^3')
    print(f'{vol_max_first = :.6g} Ang^3')
    num_alats = 12
    digits = 3
    approx_eq_alat = round(approx_eq_config.cell.cellpar()[0], digits)
    print(f'{approx_eq_alat = :.6g} Ang')

    alat_min = round(get_alat_from_V(vol_min_first, unitcell_structure), digits)
    alat_delta = round((get_alat_from_V(vol_max_first, unitcell_structure) -
        alat_min) / (num_alats - 1), digits)
    alat_max = alat_min + alat_delta * (num_alats - 1)
    alat_max_round = round(alat_max, 3)
    if abs(alat_max_round - alat_max) > 1e-6:
        # Rounding is done to to avoid not rounded numbers such as 15.088999,,,.
        print('Error: the round of alat_max is not good.')
        sys.exit()
    alat_max = alat_max_round
    vol_min = get_V_from_alat(alat_min, unitcell_structure)
    vol_max = get_V_from_alat(alat_max, unitcell_structure)

    print(f'{alat_delta = } Ang')
    print(f'{alat_min = } Ang')
    print(f'{alat_max = } Ang')
    print(f'{vol_min = :.6g} Ang^3')
    print(f'{vol_max = :.6g} Ang^3')
    alats = np.linspace(start=alat_min, stop=alat_max,
                        num=num_alats)
    alats = np.round(alats, decimals=10)  # Avoid repeating 9's at the end.
    print(f'{alats =}')
    return alats, approx_eq_config, approx_eq_alat


def make_input_for_vasp():
    """This function makes input files used by vasp for the EOS calculation."""

    # for debugging
    verbose = True
    # verbose = False

    # print(atoms_and_forces[0]['atoms'].get_positions())
    if os.path.exists('jobList_many_vol'):
        os.remove('jobList_many_vol')
    if not os.path.exists(calc_folder):
        os.mkdir(calc_folder)
    os.chdir(calc_folder)

    if not os.path.exists('b_many_vols'):
        os.mkdir('b_many_vols')
    os.chdir('b_many_vols')

    alats, __, __2 = get_alat_range(template_folder='../../template')
    # __{2}: unused variables

    for i_sample in range(num_samples):  # i_sample: index of sample
        sample_folder = get_sample_folder_name(i_sample)
        print(f'\n\n{sample_folder}')
        if not os.path.exists(sample_folder):
            os.mkdir(sample_folder)
        os.chdir(sample_folder)

        for alat in alats:
            # Pre-relaxed snapshot is used. Otherwise, many local minima
            # can be found in the energy-volume curves.
            snapshot = read(f'../../a_given_vol_relax/sample_{i_sample}/CONTCAR')
            # i_; index
            # struct_str = '{:03d}'.format(alat)    # string, padded with 0
            struct_str = f'{alat}Ang'    # string, padded with 0
            # snapshot = approx_eq_config.copy()
            snapshot.set_cell([[alat, 0, 0],
                               [0, alat, 0],
                               [0, 0, alat]], scale_atoms=True)
            folder = struct_str
            print()
            print('directory = ', folder)

            if verbose:
                print('snapshot = ', snapshot)

            if not os.path.exists(folder):
                os.mkdir(folder)
            else:
                print(f'Error: folder {folder} already exists!')
                sys.exit()
            os.chdir(folder)

            root_folder = '../../../../'
            template_folder = f'{root_folder}/template'
            write('POSCAR', snapshot, direct=True, vasp5=True)
            shutil.copy(f'{template_folder}/INCAR', '.')
            shutil.copy(f'{template_folder}/INCAR.preconverge.change', '.')
            shutil.copy(f'{template_folder}/KPOINTS', '.')
            check_link_symlink(f'{template_folder}/POTCAR', 'POTCAR')
            if use_vdw_kernel_file:
                check_link_symlink(f'{template_folder}/vdw_kernel.bindat',
                           'vdw_kernel.bindat')
            subprocess.run(f'pwd >> {root_folder}/jobList_many_vol', shell=True,
                    check=True)
            os.chdir('..')
        os.chdir('..')

    os.chdir('..')


if __name__ == '__main__':
    make_input_for_vasp()
