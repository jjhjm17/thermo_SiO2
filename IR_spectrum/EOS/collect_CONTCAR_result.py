#!/usr/bin/env python3
"""This script collects CONTCARs from equilibrium relaxation and put them as
POSCAR in the result folder."""

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
from .prepare_for_given_vol_relax import get_sample_folder_name


def collect_CONTCAR_result():
    """This function collects CONTCAR and put them in the result folder."""

    # for debugging
    # verbose = True
    verbose = False
    os.chdir(calc_folder)

    eq_vol_relax_folder = 'c_eq_vol'
    os.chdir(eq_vol_relax_folder)

    for i_sample in range(num_samples):  # i_sample: index of sample
        sample_folder = get_sample_folder_name(i_sample)
        print(f'\n\n{sample_folder}')
        os.chdir(sample_folder)

        root_folder = '../../..'

        shutil.copy('CONTCAR',
                    f'{root_folder}/{calc_folder}_result/{sample_folder}/POSCAR')
        os.chdir('..')

    os.chdir('..')
    print('CONTCARs relaxed at the equilibrium volumes are copied to the '
          'result folder as POSCARs.')


if __name__ == '__main__':
    collect_CONTCAR_result()
