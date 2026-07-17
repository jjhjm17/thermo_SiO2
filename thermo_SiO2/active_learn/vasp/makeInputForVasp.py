#!/usr/bin/env python3
"""This script makes folder for vasp input."""

import os
import sys
import shutil
import subprocess
from pathlib import Path
from collections import Counter
import numpy as np
import ase
from ase.io import read
# from mlippy.atms import ase_loadcfgs
from thermo_SiO2.mlip.read import read_SiO2
from ...mlip.read_mlip_cfg import read_mlip_cfg, read_mlip_cfg_mlippy
from ...mlip.read_config import sort_config_by_POTCAR_order
# from ...mlip.read_mlip_cfg import set_atom_symbol
# from ...util.SiO2_parameter import Si_O_H_Al_atom_symbol_tuple_mlip, Si_O_Al_atom_symbol_tuple_mlip
from ...util.SiO2_parameter import Atom_order, POTCAR_setup
from ...util.util import shell, read_output, fill_blanks
from parameters import calc_folder, atom_symbols_in_output_cfg
import parameters as param


def unique_ordered_list(seq):
    """This function returns the unique ordered list of a sequence seq."""
    # https://stackoverflow.com/questions/480214/how-do-i-remove-duplicates-from-a-list-while-preserving-order
    seen = set()
    seen_add = seen.add
    # see.add(x) is always False, but calling add(x) adds x to the set.
    return [x for x in seq if not (x in seen or seen_add(x))]


def make_POTCAR(config):
    """This function makes POTCAR from POSCAR. We assume no _sv or _pv
    POTCARS, but use the one under the symbol folder, for example,
    $pbepot/Si/POTCAR for Si."""
    all_symbols = config.get_chemical_symbols()
    unique_symbols = unique_ordered_list(all_symbols)
    # print(f'{unique_symbols =}')
    POTCAR_files = [f'$pbepot/{POTCAR_setup[symbol]}/POTCAR' for symbol in unique_symbols]
    command = 'cat ' + ' '.join(POTCAR_files) + ' > POTCAR'
    # print(f'{command =}')
    shell(command)
    return all_symbols, unique_symbols


def sort_by_chemical_formula(atoms_and_forces):
    """This function sorts atoms_and_forces object by the chemical formula and
    returns it."""
    # print(f'{atoms_and_forces[0] = }')
    formulas = [element['atoms'].get_chemical_formula() for element in
                atoms_and_forces]
    # print(f'{formulas[:3]  = }')
    index_sort = np.argsort(formulas)
    # https://numpy.org/doc/stable/reference/generated/numpy.argsort.html
    # print(f'{index_sort = }')
    sorted_atoms_and_forces = []
    for index in index_sort:
        sorted_atoms_and_forces.append( atoms_and_forces[ index])
    print('\n Configurations are sorted by the chemical formula.')

    atoms_and_forces = sorted_atoms_and_forces
    ordered_formulas = [element['atoms'].get_chemical_formula() for
                        element in atoms_and_forces]
    # print(f'{ordered_formulas = }')
    print(f'Formulas: {Counter(ordered_formulas).keys()}')
    print(f'Occurrences: {Counter(ordered_formulas).values()}\n')
    # Ref.: Count occurrences in a list
    # https://stackoverflow.com/questions/12282232/how-do-i-count-occurrence-of-unique-values-inside-a-list
    return atoms_and_forces


def get_atoms_and_forces(param):
    atom_symbol_tuple = Atom_order(atom_symbols_in_output_cfg).atom_symbol_tuple_mlip()
    if hasattr(param, 'mlip_cfg_file'):
        # atoms_and_forces = read_mlip_cfg(param.mlip_cfg_file,
        #         atom_symbol_tuple=atom_symbol_tuple,
        #         sort_method='POTCAR_order')
        # cfgs = read_mlip_cfg_mlippy(param.mlip_cfg_file,
        #         atom_symbol_tuple=param.atom_symbol_tuple)
        configs = read_mlip_cfg_mlippy(param.mlip_cfg_file,
                atom_symbol_tuple=atom_symbol_tuple)
        # atoms_and_forces = [{'atoms': cfg, 'forces':
        #                      np.zeros((len(cfg),3))} for cfg in cfgs]
    elif hasattr(param, 'POSCARs'):
        if hasattr(param, 'num_seeds'):
            if not (len(param.POSCARs) == 1 or len(param.POSCARs) == param.num_seeds):
                print("Error: The length of the variable 'POSCARs' should be 1 or num_seeds.")
                sys.exit()
            if len(param.POSCARs) == 1:
                POSCAR_files = param.POSCARs * param.num_seeds
            else:
                POSCAR_files = param.POSCARs
        else:
            POSCAR_files = param.POSCARs
        configs = []
        for POSCAR_file in POSCAR_files:
            configs.append(read(POSCAR_file, format='vasp'))
        # configs = [sort_config_by_POTCAR_order(config,
        #                                        symbols=atom_symbols_in_output_cfg) for
        #            config in configs]
        # atoms_and_forces = [{'atoms': config, 'forces':
        #                      np.zeros((len(config),3))} for config in configs]
    elif hasattr(param, 'CONTCARs'):
        CONTCAR_files = param.CONTCARs
        configs = []
        for CONTCAR_file in CONTCAR_files:
            configs.append(read(CONTCAR_file, format='vasp'))
    elif hasattr(param, 'xyz_s'):
        configs = []
        for xyz_file in param.xyz_s:
            configs.append(read(xyz_file))
        # configs = [sort_config_by_POTCAR_order(config,
        #                                        symbols=atom_symbols_in_output_cfg) for
        #            config in configs]
        # atoms_and_forces = [{'atoms': config, 'forces':
        #                      np.zeros((len(config),3))} for config in configs]
    elif hasattr(param, 'dump_s'):  # lammps dump file
        configs = []
        if hasattr(param, 'atom_symbols_input_lmp'):
            atom_symbols = param.atom_symbols_input_lmp
        else:
            atom_symbols = param.atom_symbols_in_output_cfg
        if hasattr(param, 'dump_index'):
            dump_index = param.dump_index
        else:
            dump_index = ':'
        for dump_file in param.dump_s:
            configs += read_SiO2(dump_file,
                                 atom_symbols=atom_symbols,
                                 index=dump_index)
    elif hasattr(param, 'dataf_s'):  # lammps data file
        configs = []
        if hasattr(param, 'atom_symbols_input_lmp'):
            atom_symbols = param.atom_symbols_input_lmp
        else:
            atom_symbols = param.atom_symbols_in_output_cfg
        for dataf_file in param.dataf_s:
            configs += read_SiO2(dataf_file,
                                 atom_symbols=atom_symbols,
                                 index=':')

    configs = [sort_config_by_POTCAR_order(config,
                                            symbols=atom_symbols_in_output_cfg) for
                config in configs]
    atoms_and_forces = [{'atoms': config, 'forces':
                            np.zeros((len(config),3))} for config in configs]
    if (hasattr(param, 'sort_by_chemical_formula') and
        param.sort_by_chemical_formula):
        atoms_and_forces = sort_by_chemical_formula(atoms_and_forces)
    return atoms_and_forces


def makeInputForVasp():
    # overwrite directories and files
    overwrite = False  # use only overwrite = False

    # for debugging
    verbose = True  
    # verbose = False 
    # atoms = ase_loadcfgs(mlip_cfg_file)
    # set_atom_symbol(atoms, Si_O_H_Al_atom_symbol_tuple_mlip)
    # if hasattr(param, 'atom_symbols_in_cfg') and param.atom_symbols_in_cfg == 'Si O Al':
    # if atom_symbols_in_POTCAR == 'Si O Al':
    # if atom_symbols_in_output_cfg == 'Si O Al':
    #     atom_symbol_tuple = Si_O_Al_atom_symbol_tuple_mlip
    # # elif atom_symbols_in_POTCAR == 'Si O H Al':
    # elif atom_symbols_in_output_cfg == 'Si O H Al':
    #     atom_symbol_tuple = Si_O_H_Al_atom_symbol_tuple_mlip
    # else:
    #     print('Error: please code more.')
    #     sys.exit()

    atoms_and_forces = get_atoms_and_forces(param)

    # print(atoms_and_forces[0]['atoms'].get_positions())
    if os.path.exists('jobList'):
        os.remove('jobList')
    start_dir = Path.cwd()
    if not os.path.exists(calc_folder):
        os.mkdir(calc_folder)
    os.chdir(calc_folder)

    for i_structure, this_atoms_and_forces in enumerate(atoms_and_forces):
        if hasattr(param, 'start_config_number'):
            if (i_structure < param.start_config_number) or (
                    i_structure > param.end_config_number):
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
                print(f'ERROR: Directory {folder} already exists.')
                # print('ERROR: Directory {} already exists. If you want to'.format(folder))
                # print('overwrite it, please set overwrite = True.')
                sys.exit(1)

        mkdir_if_overwrite(folder,overwrite)

        os.chdir(folder)
        shutil.copy('../../template/INCAR', '.')


        # POSCAR
        if hasattr(param, 'CONTCARs'):
            path = Path(param.CONTCARs[i_structure])
            if not path.is_absolute():
                path = start_dir / path
            shutil.copy(path, 'POSCAR')
        else:
            ase.io.write('POSCAR', snapshot, vasp5=True)

        # POTCAR
        pbepot_shell = read_output('echo $pbepot').split()
        if i_structure == 0:
            print(f'{pbepot_shell = }')
        if pbepot_shell == []:
            print('Error: please set a bash environment variable $pbepot '
                  'for the PBE POTCAR folder.')
            sys.exit()
        all_symbols, unique_symbols = make_POTCAR(snapshot)

        num_atoms = len(snapshot)
        if hasattr(param, 'nbands_less') and param.nbands_less:
            zval_s = read_output("grep ZVAL POTCAR | awk '{print $6}' ")
            zval_s = [int(float(s)) for s in zval_s.split()]
            # ZVAL_per_formula_unit = 16  # (4 + 6 * 2),  Si O2

            # Count atoms and multiply by ZVAL
            atom_counts = [all_symbols.count(sym) for sym in unique_symbols]
            nelect = np.dot(atom_counts, zval_s)
            nbands = nelect / 2 + num_atoms / 4
            # vasp manual for NBANDS
            # "In some cases, it is also possible to decrease it to NELECT/2+NIONS/4"
            
            nbands = int(np.ceil(nbands))
            fill_blanks('INCAR', blanks=['xxx__NBANDS__xxx'], variables=[f'{nbands}  # NELECT / 2 + NIONS / 4'])

        if (hasattr(param, 'preconverge') and param.preconverge):
            if os.path.isfile('../../template/INCAR.preconverge.change'):
                shutil.copy('../../template/INCAR.preconverge.change', '.')
            elif os.path.isfile('../../template/INCAR.preconverge'):
                shutil.copy('../../template/INCAR.preconverge', '.')
                fill_blanks('INCAR.preconverge', blanks=['xxx__NBANDS__xxx'],
                            variables=[f'{nbands}  # NELECT / 2 + NIONS / 4'])
            else:
                print('Error: INCAR.preconverge.change or INCAR.preconverge '
                      'are not found.')
                sys.exit()

        # os.symlink('../../template/POTCAR', 'POTCAR')

        if os.path.isfile('../../template/KPOINTS'):
            shutil.copy('../../template/KPOINTS', '.')
            # If there is no KPOINTS, then KSPACING tag is used.

        if os.path.isfile('../../template/vdw_kernel.bindat'):
            os.symlink('../../template/vdw_kernel.bindat', 'vdw_kernel.bindat')
        subprocess.run('pwd >> ../../jobList', shell=True, check=True)
        os.chdir('..')
    os.chdir('..')

if __name__ == '__main__':
    makeInputForVasp()
