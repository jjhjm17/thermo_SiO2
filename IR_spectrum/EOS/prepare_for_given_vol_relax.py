#!/usr/bin/env python3
"""This script makes folders for EOS fitting, to run vasp."""

import os
import sys
import shutil
import subprocess
import numpy as np
from ase.io import read, write
from a_parameters import (calculation_folder as calc_folder,
                          unitcell_structure, use_vdw_kernel_file, num_samples,
                          min_vol_percent, max_vol_percent,
                          POSCAR_prefix)
from ...util.util import shell


def get_sample_folder_name(index_sample):
    """This function returns the name of sample folder."""
    return f'sample_{index_sample}'


def make_input_for_vasp():
    """This function makes input files used by vasp for the EOS calculation."""

    # for debugging
    verbose = True
    # verbose = False

    # print(atoms_and_forces[0]['atoms'].get_positions())
    if os.path.exists('jobList'):
        os.remove('jobList')
    if not os.path.exists(calc_folder):
        os.mkdir(calc_folder)
    os.chdir(calc_folder)

    given_vol_relax_folder = 'a_given_vol_relax'
    if not os.path.exists(given_vol_relax_folder):
        os.mkdir(given_vol_relax_folder)
    os.chdir(given_vol_relax_folder)


    for i_sample in range(num_samples):  # i_sample: index of sample
        sample_folder = get_sample_folder_name(i_sample)
        print(f'\n\n{sample_folder}')
        if not os.path.exists(sample_folder):
            os.mkdir(sample_folder)
        os.chdir(sample_folder)

        root_folder = '../../..'
        template_folder = f'{root_folder}/template'
        snapshot = read(f'{template_folder}/{POSCAR_prefix}_{i_sample}')

        if verbose:
            print('snapshot = ', snapshot)

        write('POSCAR', snapshot, direct=True, vasp5=True)
        shell(f'cat {template_folder}/INCAR.preconverge.change '
              f' {template_folder}/INCAR > INCAR')
        shutil.copy(f'{template_folder}/KPOINTS', '.')
        os.symlink(f'{template_folder}/POTCAR', 'POTCAR')
        if use_vdw_kernel_file:
            os.symlink(f'{template_folder}/vdw_kernel.bindat',
                       'vdw_kernel.bindat')
        subprocess.run(f'pwd >> {root_folder}/jobList', shell=True, check=True)
        os.chdir('..')

    os.chdir('..')


if __name__ == '__main__':
    make_input_for_vasp()
