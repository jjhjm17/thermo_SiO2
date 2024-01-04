"""This file contains useful utility functions."""
import os
import sys
import shutil
import subprocess as sub
from datetime import datetime

def shell(command):
    """This function runs a shell command."""
    sub.run(command, shell=True, check=True)


def read_output(command):
    """Example of command: 'source ./a.defineVariables.sh; echo $scAndAtoms'"""
    return sub.run(command, shell=True, check=True,
                 capture_output=True).stdout.decode('ascii').replace('\n', ' ')


def if_exists_delete_make(folder):
    """This function deletes a folder if it exists, and makes the folder."""
    if os.path.exists(folder):
        print(f'Existing folder {folder} is deleted.')
        shutil.rmtree(folder)
    os.mkdir(folder)


def if_exists_delete_file(file, VERBOSE=True):
    """This function deletes a file if it exists."""
    if os.path.exists(file):
        if VERBOSE:
            print(f'Existing file {file} is deleted.')
        os.remove(file)


def get_lammps_random_seed():
    """This function returns a random seed for a lammps calculation."""
    # microsecond * 100
    return int(datetime.now().timestamp() * 1e6 * 100) % 899999999


def set_actual_atom_symbols(atoms, atom_symbol_tuple):
    """This function sets actual atom symbols, (0), 1, 2, 3, .. to
    atoms in atom_symbol_tuple."""
    atom_types = atoms.get_atomic_numbers()
    # print(f'{atom_types=}')
    for index in range(len(atoms)):
        atom_type = atom_types[index]  # 1, 2, ... to Si, O, ..
        atoms[index].symbol = atom_symbol_tuple[atom_type]
    return atoms


def check_link_symlink(src, dst):
    """This function checks if a symlink is good (not broken). Then make
    the symlink. src: source, dst: destination"""
    if not os.path.exists(os.path.realpath(src)):
        # realpath: Follow symlink recursively.
        print(f'\n Error: broken symlink at {src}')
        sys.exit()
    else:
        os.symlink(src, dst)
