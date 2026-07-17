"""This script reads an mlip cfg file without D4 and adds vdW D4 correction and writes the cfg file."""

import os
import sys
import shutil
import json
import subprocess
from ase.io import write
from dftd4.ase import DFTD4
from ...util.util import read_output, shell
from ...mlip.read import read_SiO2
from ...mlip.write_config import write_cfg_SiO2
from parameters import (mlip_cfg_file, output_file, functional_before,
                        atom_symbols)
import parameters
import parameters as param

# verbose = True
verbose = False

def check_computing_node():
    """This function checks if it is in computing node."""
    hostname = read_output('hostname')
    if hostname[:4] == param.server[:4].lower():
        print('Error: Are you connected to the login node? '
              'Since D4 vdW calculation is heavy, please connect to a computing node.')
        sys.exit()


# def get_D4_properties(atoms, functional):
#     """This function calculations DFT-D4 vdW using dftd4 program and returns
#     potential_energy, forces, and stresses."""
#     if os.path.exists('tmp'):
#         print(f'Error: temporary folder ./tmp exists.')
#         sys.exit()
#     os.makedirs('./tmp')
#     os.chdir('tmp')
#     write('POSCAR', atoms, format='vasp', vasp5=True)
#     # try:
#     shell(f'{param.dftd4} run -f {functional} -g grad.txt --json dftd4.json POSCAR')
#     # except subprocess.CalledProcessError:
#     #     print(f'\nError: Did you load modules using {param.dftd4_load_module}?\n')
#     with open('dftd4.json') as fin:
#         result = json.loads(fin)
#         # breakpoint()
# 
#     os.chdir('..')
#     shutil.rmtree('tmp')
#              
#     return


def add_D4():
    check_computing_node()
    if os.path.exists(output_file):
        print(f'Error: output file {output_file} exists.')
        sys.exit()

    configs = read_SiO2(mlip_cfg_file, atom_symbols)
    configs_vdW = []
    # print(f'{configs = }')
    # atoms_read = read_SiO2('tmp.cfg', atom_symbols)
    # print(f'{atoms_read = }')
    print(f'vdW energy of DFT-D4 method for {functional_before} is added.')
    # print(f'Please run manually {param.dftd4_load_module} if it is not already done.')
    print(f'Please set environments for dftd4 according to README if it is not already done.')
    print()
    print(f'Index  atoms')
    for i_structure, config in enumerate(configs):
        if hasattr(param, 'start_config_number'):
            if (i_structure < param.start_config_number) or (
                    i_structure > param.end_config_number):
                continue

        print(f'{i_structure}  {config}')
        if verbose:
            print()
            print(f'{config = }')
            print(f'{config.energy = }')
            print(f'{config.forces[:3] = }')
            print(f'{config.stresses = }')

        atoms = config
        atoms.calc = DFTD4(method=functional_before)
        disp_pot_energy = atoms.get_potential_energy()  # disp: dispersion
        disp_forces = atoms.get_forces()
        disp_stress = atoms.get_stress()
        # breakpoint()
        # get_D4_properties(config, functional_before)
        if verbose:
            print(f'{disp_pot_energy = } eV/cell')
            print(f'{disp_forces[:3] = } eV/Ang')
            print(f'{disp_stress = } eV/Ang^3')
            print(f'{atoms.get_volume() = } Ang^3')
            print(f'{disp_stress * atoms.get_volume() = } eV')

        config.energy += disp_pot_energy
        config.forces += disp_forces
        config.stresses -= disp_stress * atoms.get_volume()
        # mlip stress is - STRESS as in the vasp OUTCAR.
        configs_vdW.append(config)
        if verbose:
            print('After adding DFT-D4 values, ')
            print(f'{config.energy = }')

    write_cfg_SiO2(output_file, configs_vdW, atom_symbols=atom_symbols)
    print('File saved: ' + output_file)


    

if __name__ == '__main__':
    add_D4()
 
