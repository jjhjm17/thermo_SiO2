#!/usr/bin/env python3
"""This script creates folders and prepares jobs for H/T integration of bcc vacancy."""

import os
from glob import glob
from a_parameters import num_seeds, calc_folder
from ...util.util import if_exists_delete_file


def clean_files():
    """This function cleans calculation files."""

    os.chdir(calc_folder)
    print('Files are deleted: jobID_*, ')
    structure_file = 'cristo_192.dataf'
    links = ['SiO2_amorphous_lda_level22.mtp',
            'mlip.ini', structure_file]
    print(links)
    for index_seed in range(1, num_seeds+1):  # 1 to num_seeds
        calc_subfolder = f'seed_{str(index_seed).zfill(3)}'
        print(calc_subfolder)
        os.chdir(calc_subfolder)


        jobID_files = glob('jobID_*')
        for jobID_file in jobID_files:
            if_exists_delete_file(jobID_file, VERBOSE=False)
        for link in links:
            if_exists_delete_file(link, VERBOSE=False)
        os.chdir('..')
    os.chdir('..')


if __name__ == "__main__":
    clean_files()
