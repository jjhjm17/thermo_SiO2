#!/usr/bin/env python3
"""This script obtains energy volume curve from the calculated DFT result."""

import os
import sys
import numpy as np
# from python_fit.fitToEOS import fitToEOS
from a_parameters import (calculation_folder as calc_folder,
                          unitcell_structure, num_samples)
from ...util.util import read_output
from .make_folder import get_alat_range, get_sample_folder_name


def get_E_V_curve():
    """This function obtains E-V curve."""

    if os.path.exists('E_V_input'):
        print('The existing E_V_input file is removed.')
        os.remove('E_V_input')
    os.chdir(calc_folder)

    alats, approx_eq_config = get_alat_range()
    num_atoms = len(approx_eq_config)

    for i_sample in range(num_samples):  # i_sample: index of sample
        V_per_atom = []
        E_per_atom = []
        sample_folder = get_sample_folder_name(i_sample)
        print(f'\n\n{sample_folder}')
        os.chdir(sample_folder)

        for alat in alats:
            # i_; index
            # struct_str = '{:03d}'.format(alat)    # string, padded with 0
            struct_str = f'{alat}Ang'    # string, padded with 0
            folder = struct_str
            # print('directory = ', folder)
            os.chdir(folder)

            if unitcell_structure != 'amorphous_cubic':
                print("Error: for now only unitcell_structure = 'amorphous_cubic' is"
                      " supported.")
                sys.exit()

            V_per_atom.append(alat**3 / num_atoms)  # Ang^3/atom
            E_per_atom.append( float(read_output('tail -n1 OSZICAR').split()[4])
                              / num_atoms)

            os.chdir('..')
        np.savetxt('E_V_input', np.column_stack((V_per_atom, E_per_atom)),
                   fmt='%.16g', header='volume (Ang^3/atom)  energy (eV/atom)')
        os.chdir('..')
    os.chdir('..')

    # fitToEOS(volumeUnit='Ang^3', energyUnit='eV', file='E_V_input',
    #          latType='sc', cBya=1, fitType='Vinet',  # sc: simple cubic
    #         $Vdef*1., $Bdef*1., $Bmin*1.,
    #     $Bmax*1., $BderDef*1., $BderMin*1.,
    #     $BderMax*1., $mesh)

    print('File written: E_V_input')



if __name__ == '__main__':
    get_E_V_curve()
