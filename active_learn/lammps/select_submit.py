#!/usr/bin/env python3
"""This script collects preselected configurations and run a selection job."""
import os
import subprocess
import shutil
from glob import glob
from ...util.util import (shell, if_exists_delete_file, if_exists_delete_make,
                          fill_blanks)
from a_parameters import calc_folder, num_seeds, server
import a_parameters
import a_parameters as param

def collect_preselected():
    """This function collects preselected configurations."""
    if_exists_delete_file('all_sampled.cfg')
    if hasattr(a_parameters, 'select_start_seed'):
        select_start_seed = a_parameters.select_start_seed
        select_end_seed = a_parameters.select_end_seed
        print(f'{select_start_seed = }, {select_end_seed = }')
        shell(f'touch all_sampled.cfg')
        for index_seed in range(num_seeds):
            if (index_seed < select_start_seed or select_end_seed <
                index_seed):
                continue
            calc_subfolder = f'{str(index_seed).zfill(3)}'
            # preselected_files = f'{calc_folder}/{calc_subfolder}/preselected.cfg'
            preselected_files = sorted(glob(f'{calc_folder}/{calc_subfolder}/preselected*.cfg'))
            # if os.path.isfile(preselected_file):
            if len(preselected_files) >= 1:
                # if calculation finished good and wrote the preselected file
                for pre_file in preselected_files:
                    # pre_file: preselected file
                    shell(f'cat {pre_file}  >> all_sampled.cfg')
            else:
                print(f'Warning: {calc_subfolder} does not have '
                      'preselected*.cfg file.')
    else:
        shell(f'cat {calc_folder}/*/preselected*.cfg > all_sampled.cfg')

    print('File written: all_sampled.cfg')
    print('Number of configurations: ', end='')
    shell(f'grep END_CFG all_sampled.cfg | wc --lines')


def select_submit():
    """This function runs a selection job."""
    collect_preselected()
    if_exists_delete_make('select')
    os.chdir('select')
    # shutil.copy('../template/select.sh', '.')
    print(f'{server = }')
    shutil.copy(f'../template/select.{server}.sh', '.')

    if hasattr(a_parameters, 'almtp'):
        fill_blanks(file=f'select.{server}.sh',
                    blanks=['xxx__almtp__xxx', 'xxx__train_cfg__xxx'],
                    variables=[a_parameters.almtp, a_parameters.train_cfg])
    shell(f'sbatch select.{server}.sh')
    print('In folder select, job is submitted.')
    os.chdir('..')

if __name__ == "__main__":
    select_submit()
