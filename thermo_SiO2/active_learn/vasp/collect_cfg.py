#!/usr/bin/env python3
"""This script collects calculated configuration to cfg file."""

import os
import shutil
import sys
from time import sleep
# import numpy as np
from ase.io import read
from parameters import (calc_folder, output_file, mlip_version,
                        server, atom_symbols_in_output_cfg)
import parameters as param
from ...util.util import shell, read_output
from ...mlip.read_mlip_cfg import read_mlip_cfg
from ...mlip.read_mlip_cfg import read_mlip_cfg_mlippy
# from mlippy.atms import loadcfgs, savecfgs
# from ...mlip.read_mlip_cfg import set_atom_symbol
from ...util.SiO2_parameter import Si_O_H_Al_atom_symbol_tuple_mlip
if hasattr(param, 'added_train_cfg'):
    sys.path.insert(0, '../a.lammps')
    import a_parameters as parameters_lammps
from thermo_SiO2.active_learn.vasp.makeInputForVasp import get_atoms_and_forces


def collect_cfg():
    """This function collects calculated configurations."""

    if os.path.exists(output_file):
        print(f'Error: output file {output_file} exists.')
        sys.exit()

    # atoms = ase_loadcfgs(mlip_cfg_file)
    # set_atom_symbol(atoms, Si_O_H_Al_atom_symbol_tuple_mlip)
    if hasattr(param, 'mlip_cfg_file'):
        # atoms_and_forces = read_mlip_cfg(param.mlip_cfg_file,
        #         atom_symbol_tuple=Si_O_H_Al_atom_symbol_tuple_mlip,
        #         sort_method='POTCAR_order')
        cfgs = read_mlip_cfg_mlippy(param.mlip_cfg_file,
                atom_symbol_tuple=Si_O_H_Al_atom_symbol_tuple_mlip)
        # here, Si_O_H_Al_atom_symbol_tuple_mlip is not used later.

    # atoms_and_forces = get_atoms_and_forces(param)

    os.chdir(calc_folder)


    # for i_structure in range( len( atoms_and_forces)):
    for i_structure in range( len( cfgs)):
        if hasattr(param, 'start_config_number'):
            if (i_structure < param.start_config_number) or (
                    i_structure > param.end_config_number):
                continue
        # i_; index
        # struct_str = '{:03d}'.format(i_structure)    # string, padded with 0
        struct_str = f'{i_structure:04d}'    # string, padded with 0
        folder = struct_str

        if i_structure == 0:
            # ISMEAR = int(read_output(f"grep 'ISMEAR =' {folder}/OUTCAR").split(sep=';')[0].split(sep='=')[1])
            ISMEAR = int(read_output(f"grep 'ISMEAR =    ' {folder}/OUTCAR").split(sep=';')[0].split(sep='=')[1])
            if ISMEAR == 0:  # Gaussian
                # 'mlp_par_fits_to_energy_based_on_ismear' is usually 
                # used for reading energy. For MP smearing,
                # read E0. For the Gaussian smearing, read F.
                # For the Gaussian smearing, we use mlp for speed.
                if server == 'justus':
                    mlp_binary = 'mlp'
                elif server == 'fritz':
                    mlp_binary = 'mlp2ser'
            else:
                mlp_binary = 'mlp_par_fits_to_energy_based_on_ismear'
            print(f'{ISMEAR = }, {mlp_binary = }')
            sleep(5)  # sec

        print(f'directory = {folder}   ', end='')

        os.chdir(folder)


        if mlip_version == 2:
            print('Error: mlip_version = 2 is not supported for remapping atom types.')
            sys.exit()
            # shell(f'{mlp_binary} '
            #     f'convert-cfg OUTCAR ../../{output_file} --append '
            #     '--input-format=vasp-outcar')
        elif mlip_version == 3:
            if server == 'fritz':
                mlp_binary = 'mlp3ser'
            elif server == 'justus':
                mlp_binary = 'mlp3'
            else:
                print('Error: please code more for this server.')
                sys.exit()

            shell(f'{mlp_binary} '
                f'convert OUTCAR before_remap.cfg '
                '--input_format=outcar')
            # with open('POSCAR', 'r') as fin:
            #     for index, line in enumerate(fin):
            #         if index == 5:
            #             atom_symbols_in_POTCAR = line.split()  #replace('  ', ' ').strip()
            #         elif index > 5:
            #             break
            read_outcar = read('OUTCAR')
            atom_symbols_in_POTCAR = list(read_outcar.symbols.indices().keys())
            # print(f'{atom_symbols_in_POTCAR = }')
            symbols_out = atom_symbols_in_output_cfg.split()
            remap_numbers = [f'{symbols_out.index(symbol)}' for
                             symbol in atom_symbols_in_POTCAR]
            # If remap_numbers is too short, mlp3 does not work. So we append
            # more numbers without meaning to keep the length as expected.
            if len(remap_numbers) < len(symbols_out):
                full_list = [f'{x}' for x in list(range(len(symbols_out)))]
                remaining = [x for x in full_list if x not in remap_numbers]
                remap_numbers.extend(remaining)
            remap_numbers = ' '.join(remap_numbers)
            # print(f'{remap_numbers = }')
            shell(f'{mlp_binary} '
                f'remap_species before_remap.cfg remapped.cfg {remap_numbers}')
            os.remove('before_remap.cfg')
            shell(f'cat remapped.cfg >> ../../{output_file}')
            os.remove('remapped.cfg')

            # if atom_symbols_in_output_cfg == 'Si O H Al':
            #     if atom_symbols_in_POTCAR == 'O H':
            #         shell(f'{mlp_binary} '
            #             f'remap_species before_remap.cfg remapped.cfg 1 2 ')
            #         # 1 : O, 2 : H
            #         os.remove('before_remap.cfg')
            #     elif atom_symbols_in_POTCAR == 'Si O Al':
            #         shell(f'{mlp_binary} '
            #             f'remap_species before_remap.cfg remapped.cfg 0 1 3 ')
            #         # 0: Si, 1 : O, 3 : Al
            #         os.remove('before_remap.cfg')
            #     elif atom_symbols_in_POTCAR == atom_symbols_in_output_cfg:
            #         shutil.move('before_remap.cfg', 'remapped.cfg')
            #     else:
            #         print('Error: please code more for atom symbols (No. 2).')
            #         sys.exit()
            #     shell(f'cat remapped.cfg >> ../../{output_file}')
            #     os.remove('remapped.cfg')
            # elif atom_symbols_in_output_cfg == 'Si O Al':
            #     if atom_symbols_in_POTCAR == atom_symbols_in_output_cfg:
            #         shutil.move('before_remap.cfg', 'remapped.cfg')
            #     else:
            #         print('Error: please code more for atom symbols (No. 3).')
            #         sys.exit()
            #     shell(f'cat remapped.cfg >> ../../{output_file}')
            # else:
            #     print('Error: please code more for these atom symbols.')
            #     sys.exit()
        else:
            print('Error: mlip_version is not supported.')
            sys.exit()
        # shell(f'mlp3ser convert OUTCAR ../../{output_file} --append '
        #     '--input_format=outcar')
        os.chdir('..')
    os.chdir('..')

    atoms_and_forces_out = read_mlip_cfg(output_file,
            atom_symbol_tuple=Si_O_H_Al_atom_symbol_tuple_mlip,
            sort_method='POTCAR_order')
    print(f'Total {len(atoms_and_forces_out)} configurations are collected in '
          f'this file: {output_file}')

    if hasattr(param, 'added_train_cfg'):
        added_train_cfg = param.added_train_cfg
        shell(f'cat {parameters_lammps.train_cfg} {output_file} > {added_train_cfg}')
        print(f'Added train set is written: {added_train_cfg}, total configurations: ', end='')
        shell(f'grep END_CFG {added_train_cfg } | wc --lines')

    # # correct atom type
    # if (hasattr(parameters, 'atom_symbols_in_cfg') and
    #     parameters.atom_symbols_in_cfg == 'Si O Al'):
    #     cfgs = loadcfgs(output_file)
    #     convert_type = {0:0, 1:1, 2:3}  # Change Al type from 2 to 3.
    #     for cfg in cfgs:
    #         print(f'Before, {cfg.types = }')
    #         cfg.types = np.array([convert_type[atom_type] for atom_type in
    #                               cfg.types])
    #         print(f'After, {cfg.types = }')
    #     os.remove(output_file)
    #     savecfgs(output_file, cfgs)

if __name__ == '__main__':
    collect_cfg()
