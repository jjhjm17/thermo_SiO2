#!/usr/bin/env python3
"""This script collects preselected configurations and run a selection job."""
import os
import subprocess
import shutil
from ...util.util import shell, if_exists_delete_file, if_exists_delete_make
from a_parameters import calc_folder

def collect_preselected():
    """This function collects preselected configurations."""
    if_exists_delete_file('all_sampled.cfg')
    shell(f'cat {calc_folder}/*/preselected.cfg > all_sampled.cfg')
    print('File written: all_sampled.cfg')

def select_submit():
    """This function runs a selection job."""
    collect_preselected()
    if_exists_delete_make('select')
    os.chdir('select')
    shutil.copy('../template/select.sh', '.')
    shell('sbatch select.sh')
    print('In folder select, job is submitted.')
    os.chdir('..')

if __name__ == "__main__":
    select_submit()
