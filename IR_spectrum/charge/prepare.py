#!/usr/bin/env python3
"""This script makes folders for charge calculation, to run vasp."""

import os
import shutil
from ase.io import read, write
from a_parameters import (calculation_folder as calc_folder,
                          num_samples, use_vdw_kernel_file, symbols)
from ..EOS.util import get_sample_folder_name
from ...mlip.read_config import read_SiO2_dump, sort_config_by_POTCAR_order
from ...util.util import shell
from ...util.SiO2_parameter import Si_O_Al_atom_symbol_tuple_lammps, Si_O_H_Al_atom_symbol_tuple_lammps
import a_parameters


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
        if hasattr(a_parameters, 'POSCAR_files'):
            POSCAR_files = a_parameters.POSCAR_files
            POSCAR_file = f'{root_folder}/{POSCAR_files[i_sample]}'
            snapshot = read(POSCAR_file)
        elif hasattr(a_parameters, 'dump_files'):
            if symbols == 'Si O H Al':
                atom_symbol_tuple = Si_O_H_Al_atom_symbol_tuple_lammps
            elif symbols == 'Si O Al':
                atom_symbol_tuple = Si_O_Al_atom_symbol_tuple_lammps
            else:
                print('Error: please define the symbols variable correctly.')
                sys.exit()
            dump_file = f'{root_folder}/{a_parameters.dump_files[i_sample]}'
            snapshot = read_SiO2_dump(dump_file, index='-1',
                                      atom_symbol_tuple=atom_symbol_tuple)[0]
            if symbols == 'Si O H Al':
                atom_symbol_tuple = Si_O_H_Al_atom_symbol_tuple_lammps
                snapshot = sort_config_by_POTCAR_order(snapshot)
            else:
                print('Atoms are not sorted by the POTCAR order')




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
