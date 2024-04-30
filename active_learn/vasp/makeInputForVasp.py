#!/usr/bin/env python3
"""This script makes folder for vasp input."""

import os
import sys
import shutil
import subprocess
import ase
# from mlippy.atms import ase_loadcfgs
from ...mlip.read_mlip_cfg import read_mlip_cfg
# from ...mlip.read_mlip_cfg import set_atom_symbol
from ...util.SiO2_parameter import Si_O_H_Al_atom_symbol_tuple_mlip, Si_O_Al_atom_symbol_tuple_mlip
from parameters import mlip_cfg_file, calc_folder
import parameters


def makeInputForVasp():
    # overwrite directories and files
    overwrite = False  # use only overwrite = False

    # for debugging
    verbose = True  
    # verbose = False 
    # atoms = ase_loadcfgs(mlip_cfg_file)
    # set_atom_symbol(atoms, Si_O_H_Al_atom_symbol_tuple_mlip)
    if hasattr(parameters, 'atom_symbols_in_cfg') and parameters.atom_symbols_in_cfg == 'Si O Al':
        atom_symbol_tuple = Si_O_Al_atom_symbol_tuple_mlip
    else:
        atom_symbol_tuple = Si_O_H_Al_atom_symbol_tuple_mlip
    atoms_and_forces = read_mlip_cfg(mlip_cfg_file,
            atom_symbol_tuple=atom_symbol_tuple,
            sort_method='POTCAR_order')

    # print(atoms_and_forces[0]['atoms'].get_positions())
    if os.path.exists('jobList'):
        os.remove('jobList')
    if not os.path.exists(calc_folder):
        os.mkdir(calc_folder)
    os.chdir(calc_folder)

    for i_structure, this_atoms_and_forces in enumerate(atoms_and_forces):
        if hasattr(parameters, 'start_config_number'):
            if (i_structure < parameters.start_config_number) or (
                    i_structure > parameters.end_config_number):
                continue
        # i_; index
        # snapshot = ase.io.read('config_step_0.cfg')
        struct_str = '{:04d}'.format(i_structure)    # string, padded with 0
        snapshot = this_atoms_and_forces['atoms']

        folder = struct_str
        print()
        print('directory = ', folder)
        if verbose: 
            print('snapshot = ', snapshot)
            # print('snapshot.get_positions()[0:10] = ', snapshot.get_positions()[0:10])
            # print("'{:.15f}'.format(snapshot.get_scaled_positions()[1][1]) = ", '{:.15f}'.format(snapshot.get_scaled_positions()[1][1]))

        def mkdir_if_overwrite(folder, overwrite):
            if overwrite:
                if os.path.exists(folder):
                    shutil.rmtree(folder)
            try:
                os.mkdir(folder) 
            except FileExistsError:
                print('ERROR: Directory {} already exists.')
                # print('ERROR: Directory {} already exists. If you want to'.format(folder))
                # print('overwrite it, please set overwrite = True.')
                sys.exit(1)

        mkdir_if_overwrite(folder,overwrite)

        os.chdir(folder)
        ase.io.write('POSCAR', snapshot, vasp5=True)
        # num_atoms = len(snapshot)
        # ZVAL_per_formula_unit = 16  # (4 + 6 * 2),  Si O2
        # NELECT = num_atoms / 3 * ZVAL_per_formula_unit  # 1024 for 192 atoms
        # NBANDS = NELECT / 2 + num_atoms / 4  # See explanation in INCAR.
        # Value for SiO2.
        # The number of approximate unoccupied bands is 50 for 200 atom SiO2.
        # print(f'{num_atoms=}')
        # print(f'{NELECT =}')
        # print(f'{NBANDS =}')
        # NBANDS = int(np.ceil(NBANDS))
        # with open('../../INCAR', 'r') as fin:
        #     with open('./INCAR', 'w') as fout:
        #         for line in fin:
        #             fout.write(line.replace('xxxNBANDSxxx',
        #                 f'{NBANDS}  # NELECT = {NELECT}, NIONS = {num_atoms}'))
        shutil.copy('../../template/INCAR', '.')
        if (hasattr(parameters, 'preconverge') and parameters.preconverge ==
                True):
            shutil.copy('../../template/INCAR.preconverge.change', '.')
        shutil.copy('../../template/KPOINTS', '.')
        os.symlink('../../template/POTCAR', 'POTCAR')
        os.symlink('../../template/vdw_kernel.bindat', 'vdw_kernel.bindat')
        subprocess.run('pwd >> ../../jobList', shell=True, check=True)
        os.chdir('..')
    os.chdir('..')

if __name__ == '__main__':
    makeInputForVasp()
