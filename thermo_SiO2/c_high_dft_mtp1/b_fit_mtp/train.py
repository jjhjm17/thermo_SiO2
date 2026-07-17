#!/usr/bin/env python3
"""This script generate a folder and submit job to train the 1/10 of 
allcfgs. After the 1st training, jobs for 2nd and 3rd trainings can be
submitted.
"""

import os
import os.path as path
import subprocess as sub
from ...util.util import shell, if_exists_delete_make
from .server_parameters import init_mtp_folder, 
from d_defineMore import (min_dist, max_dist, mtp_levels, radial_basis_size,
    job_file, species_count)

for mtp_level in mtp_levels:
    print('MTP level:', mtp_level)
    if_exists_delete_make(f'mtp_level_{mtp_level}')
    os.chdir(f'mtp_level_{mtp_level}')

    # if_exists_delete_make('mlp_a_every_10')
    # os.chdir('mlp_a_every_10')

    # shell('ln -s ../training_set.cfg')
    shell(f'cp {init_mtp_folder}/{mtp_level}.mtp init.mtp')
    with open('init.mtp', 'r') as fin:
        with open('init_modified.mtp', 'w') as fout:
            for line in fin:
                line = line.replace('species_count = 1',
                    f'species_count = {species_count}')
                line = line.replace('min_dist = 2', f'min_dist = {min_dist}')
                line = line.replace('max_dist = 5', f'max_dist = {max_dist}')
                line = line.replace('radial_basis_size = 8',
                       f'radial_basis_size = {radial_basis_size}')
                fout.write(line)
    shell('mv init_modified.mtp init.mtp')
    shell(f'echo "{mtp_level}" > mtp_level')
    
    job_file_path = path.dirname(path.abspath(__file__)) + f'/{job_file}'
    shell(f'cp {job_file_path} .')
    #  shell('sbatch run.mlip_24cores.sh . | tee job_ID')
    shell(f'sbatch {job_file} . | tee job_ID')
    # shell('./run.mlip_24cores.sh')  # debug
    os.chdir('..')


