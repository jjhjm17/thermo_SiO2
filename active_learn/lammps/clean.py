#!/usr/bin/env python3
"""This script cleans duplicate files after calculation."""
import os
from glob import glob
from a_parameters import calc_folder, num_seeds

def clean():
    """This function runs a selection job."""
    os.chdir(calc_folder)
    print('Files in the calculation folder will be removed: in.file, *.almtp, jobID_*')

    for index_seed in range(num_seeds):
        calc_subfolder = f'{str(index_seed).zfill(3)}'
        os.chdir(calc_subfolder)
        try:
            os.remove('in.file')
        except FileNotFoundError:
            pass
        mtp_file = glob('*.almtp')
        if len(mtp_file) > 0:
            os.unlink(mtp_file[0])
        jobID_file = glob('jobID_*')
        if len(jobID_file) > 0:
            os.unlink( jobID_file[0])

        os.chdir('..')
    os.chdir('..')


if __name__ == "__main__":
    clean()
