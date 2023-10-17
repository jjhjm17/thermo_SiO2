#!/usr/bin/env python3

import sys
from .view_SiO2_mlip  import view_SiO2_mlip
from .view_SiO2_lammps import view_SiO2_lammps

def view_SiO2(file):
    """This function shows lammps or mlip SiO2 file."""
    # split = file.split('.')
    # if len(split) == 0:
    # suffix = file.split('.')[-1]
    # suffix_2nd = file.split('.')[-2]
    # if suffix == 'cfg' or suffix_2nd == 'cfg':
    if 'cfg' in file:
        view_SiO2_mlip(file)
    else:
        view_SiO2_lammps(file)

if __name__ == '__main__':
    view_SiO2(sys.argv[1])
