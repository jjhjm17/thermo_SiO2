#!/usr/bin/env python3
"""This script obtains BORN charge file from vasp output."""

import os
import shutil
import numpy as np
from ase.io import read
from a_parameters import (calculation_folder as calc_folder,
                          num_samples, use_vdw_kernel_file)
from ...util.util import shell
from ..EOS.util import get_sample_folder_name


def make_input_for_vasp():
    """This function makes input files used by vasp for the charge calculation."""

    # for debugging
    verbose = True
    # verbose = False

    os.chdir(calc_folder)

    for i_sample in range(num_samples):  # i_sample: index of sample
        sample_folder = get_sample_folder_name(i_sample)
        print(f'\n\n{sample_folder}')
        os.chdir(sample_folder)
        shell('phonopy-vasp-born > BORN')
        os.chdir('..')

    os.chdir('..')

    result_folder = f'{calc_folder}_result'
    if not os.path.exists(result_folder):
        os.mkdir(result_folder)
    os.chdir(result_folder)
    # charges = []
    for i_sample in range(num_samples):  # i_sample: index of sample
        sample_folder = get_sample_folder_name(i_sample)
        charge_file = f'{sample_folder}_BORN'
        if os.path.exists(charge_file):
            print('Error: charge file already exists.')
            sys.exit()
        root_folder = '..'
        shutil.move(f'{root_folder}/{calc_folder}/{sample_folder}/BORN',
                    charge_file)
        shutil.copy(f'{root_folder}/{calc_folder}/{sample_folder}/POSCAR',
                    f'{sample_folder}_POSCAR')
        # # average
        # charges.append(np.loadtxt(charge_file))
    # charge_avg = np.average(charges, axis=0)
    # charge_std = np.std(charges, axis=0)
    # with open(f'{sample_folder}_BORN', 'r') as fin:
    #     header = fin.read().splitlines()[0]
    # np.savetxt('BORN_avg', charge_avg, header=header, fmt='%.8g')
    # np.savetxt('BORN_std', charge_std, header=header, fmt='%.8g')


    # print(f'Folder made: {calc_folder}_result with BORN, POSCAR, BORN_avg.')
    print(f'Folder made: {calc_folder}_result with BORN, POSCAR.')


if __name__ == '__main__':
    make_input_for_vasp()
