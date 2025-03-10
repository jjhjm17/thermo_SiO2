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


def get_lammps_random_seed(rng=None):
    """This function returns a random seed for a lammps calculation.
    rng: numpy random number generator
    ex) rng = np.random.default_rng(12345)
    The random number generator makes reproducible numbers from a seed."""
    rand_max = 899999999
    if rng == None:
        rand = int(datetime.now().timestamp() * 1e6 * 100) % rand_max
    else:
        rand = rng.integers(low=0, high=rand_max, size=1)[0]
    # microsecond * 100
    return rand


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


def fill_blanks(file, blanks, variables):
    """This function fills blanks by variables.
    ex) file = 'in.file'
    blanks = ['xxx__seed__xxx']
    variables = [seed]"""
    if len(variables) != len(blanks):
        print('Error: the length of variables and blanks should be '
              'same.')
        sys.exit()
    with open(file, 'r') as fin:
        in_file_lines = fin.readlines()
    with open(f'{file}_new', 'w') as out_file:
        for line in in_file_lines:
            for variable, blank in zip(variables, blanks):
                line = line.replace(blank, variable)
            out_file.write(line)
    os.remove(file)
    os.rename(f'{file}_new', file)
