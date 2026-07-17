#!/usr/bin/env python3
"""This module reads OUTCARS and writes configurations."""
import sys  # debug
import os.path
import numpy as np
import a_parameters
from a_parameters import outcar_files, write_initial_training_set
from .server_parameters import mlp
from ...util.util import shell, read_output

def vaspOUTCAR_to_cfg(outcar_files):
    """This function changes vasp OUTCAR to mlip cfg file to allcfgs.cfg.
    If initial_training_every > 1, every 'initial_training_every' configurations are written to
    training_set.cfg, which is an initial training set."""
    if hasattr(a_parameters, 'select_index_file'):
        do_select_index = True
        selected_indexes = np.loadtxt(a_parameters.select_index_file,
                dtype=np.int_)
        index_of_selection = 0
    if hasattr(a_parameters, 'write_initial_training_set'):
        initial_training_every = a_parameters.initial_training_every

    print(f'{initial_training_every = }')
    if os.path.isfile('allcfgs.cfg'):
        print("The existing allcfgs.cfg file is removed.\n")
        os.remove('allcfgs.cfg')
    if os.path.isfile('training_set.cfg'):
        print("The existing training_set.cfg file is removed.\n")
        os.remove('training_set.cfg')

    for index, file in enumerate(outcar_files):
        print(f'{index}  ', end='')
        if do_select_index:
            if index == selected_indexes[index_of_selection]:
                print(f'Index {index} is selected. ', end='')
                shell(f'{mlp} convert-cfg --input-format=vasp-outcar --append {file} allcfgs.cfg')
                index_of_selection += 1
                if index_of_selection == len(selected_indexes):
                    print('All configurations are selected.')
                    break
            else:
                continue
        else:
            shell(f'{mlp} convert-cfg --input-format=vasp-outcar --append {file} allcfgs.cfg')

    num_cfgs = read_output('grep Size allcfgs.cfg | wc -l')
    # word count, count the number of lines
    print(f"\nNumber of configurations: {num_cfgs}")

    if write_initial_training_set:
        shell(f'{mlp} subsample allcfgs.cfg training_set.cfg {initial_training_every}')  #
        # every initial_training_every configs

    print()
    mindist = read_output(f'{mlp} mindist allcfgs.cfg')
    print(f"Minimum distance (Ang): {mindist}")

if __name__ == '__main__':
    vaspOUTCAR_to_cfg(outcar_files)
