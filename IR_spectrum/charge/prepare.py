#!/usr/bin/env python3
"""This script makes folders for charge calculation, to run vasp."""

import os
import shutil
from ase.io import read, write
from a_parameters import (calculation_folder as calc_folder,
                          num_samples, use_vdw_kernel_file,
                          POSCAR_files)
from ...util.util import shell
from ..EOS.util import get_sample_folder_name


def make_input_for_vasp():
    """This function makes input files used by vasp for the charge calculation."""

    # for debugging
    verbose = True
    # verbose = False

    # print(atoms_and_forces[0]['atoms'].get_positions())
    if os.path.exists('jobList'):
        os.remove('jobList')
    if not os.path.exists(calc_folder):
        os.mkdir(calc_folder)
    os.chdir(calc_folder)

    for i_sample in range(num_samples):  # i_sample: index of sample
        sample_folder = get_sample_folder_name(i_sample)
        print(f'\n\n{sample_folder}')
        if not os.path.exists(sample_folder):
            os.mkdir(sample_folder)
        os.chdir(sample_folder)

        root_folder = '../..'
        template_folder = f'{root_folder}/template'
        shutil.copy(f'{template_folder}/INCAR', '.')
        # POSCAR_file = f'{root_folder}/{POSCAR_folder}/{sample_folder}/POSCAR'
        POSCAR_file = f'{root_folder}/{POSCAR_files[i_sample]}'
        snapshot = read(POSCAR_file)

        if verbose:
            print('snapshot = ', snapshot)

        # shutil.copy(POSCAR_file, '.')
        write('POSCAR', snapshot, format='vasp')
        shutil.copy(f'{template_folder}/KPOINTS', '.')
        os.symlink(f'{template_folder}/POTCAR', 'POTCAR')
        if use_vdw_kernel_file:
            os.symlink(f'{template_folder}/vdw_kernel.bindat',
                       'vdw_kernel.bindat')
        shell(f'pwd >> {root_folder}/jobList')
        os.chdir('..')

    os.chdir('..')


if __name__ == '__main__':
    make_input_for_vasp()
