#!/usr/bin/env python3
"""This script obtains the maximum forces in configurations."""
from glob import glob
from ...mlip.get_properties import get_max_force
from parameters import output_file, atom_symbols_in_output_cfg
import parameters

if __name__ == '__main__':

    print(f'File: {output_file}')
    get_max_force(output_file, atom_symbols_in_output_cfg, 'fig-force-new.pdf')
    if hasattr(parameters, 'added_train_cfg'):
        added_train_cfg = parameters.added_train_cfg
    else:
        added_train_cfg = glob('train_*_added.cfg')[0]
    print()
    print(f'File: {added_train_cfg}')
    get_max_force(added_train_cfg, atom_symbols_in_output_cfg, 'fig-force-added.pdf')
