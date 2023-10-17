#!/usr/bin/env python3
"""This script makes folder for vasp input. NBANDS is adjusted."""

import os
import sys
import shutil
import subprocess 
import numpy as np
import ase
from ..mlip import read_mlip_cfg
from parameters import mlip_cfg_file, calc_folder

def makeInputForVasp():
    # overwrite directories and files
    # overwrite = True  # debug
    overwrite = False

    # for debugging
    verbose = True  
    # verbose = False 
    atoms_and_forces = read_mlip_cfg.read_mlip_cfg(mlip_cfg_file,
            atom_symbol_tuple={0: 'Si', 1: 'O'})  # SiO2

    # print(atoms_and_forces[0]['atoms'].get_positions())
    if os.path.exists('jobList'):
        os.remove('jobList')
    if not os.path.exists(calc_folder):
        os.mkdir(calc_folder)
    os.chdir(calc_folder)

    for i_structure, this_atoms_and_forces in enumerate(atoms_and_forces):
        # i_; index
        # snapshot = ase.io.read('config_step_0.cfg')
        struct_str = '{:05d}'.format(i_structure)    # padded with 0
        snapshot = this_atoms_and_forces['atoms']

        folder = 'vasp_struct_' + struct_str
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
                print('ERROR: Directory {} already exists. If you want to'.format(folder))
                print('overwrite it, please set ovewrite = True.')
                sys.exit(1)

        mkdir_if_overwrite(folder,overwrite)

        os.chdir(folder)
        ase.io.write('POSCAR', snapshot, vasp5=True)
        num_atoms = len(snapshot)
        ZVAL_per_formula_unit = 16  # (4 + 6 * 2),  Si O2
        NELECT = num_atoms / 3 * ZVAL_per_formula_unit  # 1024 for 192 atoms
        NBANDS = NELECT / 2 + num_atoms / 4  # See explanation in INCAR.
        # Value for SiO2.
        # The number of approximate unoccupied bands is 50 for 200 atom SiO2.
        # print(f'{num_atoms=}')
        # print(f'{NELECT =}')
        # print(f'{NBANDS =}')
        NBANDS = int(np.ceil(NBANDS))
        with open('../../INCAR', 'r') as fin:
            with open('./INCAR', 'w') as fout:
                for line in fin:
                    fout.write(line.replace('xxxNBANDSxxx',
                        f'{NBANDS}  # NELECT = {NELECT}, NIONS = {num_atoms}'))
        shutil.copy('../../INCAR.preconverge.change', '.')
        shutil.copy('../../KPOINTS', '.')
        os.symlink('../../POTCAR', 'POTCAR')
        subprocess.run('pwd >> ../../jobList', shell=True, check=True)
        os.chdir('..')
    os.chdir('..')

if __name__ == '__main__':
    makeInputForVasp()
