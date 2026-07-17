#!/usr/bin/env python3
"""This module reads cfg files and subsamples only 1 every
initial_training_every (usually 10) configurations."""
import sys  # debug
import os.path
import shutil
import numpy as np
import a_parameters
from a_parameters import input_cfg_file, server 
from .server_parameters import mlp
from ...util.util import shell, read_output

def subsample(input_cfg_file):
    """This function changes vasp OUTCAR to mlip cfg file to allcfgs.cfg.
    If initial_training_every > 1, every 'initial_training_every' configurations are written to
    training_set.cfg, which is an initial training set."""
    if hasattr(a_parameters, 'write_initial_training_set'):
        write_initial_training_set = a_parameters.write_initial_training_set
        initial_training_every = a_parameters.initial_training_every

    print(f'{initial_training_every = }')
    if os.path.isfile('allcfgs.cfg'):
        print("The existing allcfgs.cfg file is removed.\n")
        os.remove('allcfgs.cfg')
    if os.path.isfile('training_set.cfg'):
        print("The existing training_set.cfg file is removed.\n")
        os.remove('training_set.cfg')

    shutil.copy(input_cfg_file, 'allcfgs.cfg')
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
    # subsample(outcar_files)
    subsample(input_cfg_file)
