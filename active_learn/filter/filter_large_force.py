#!/usr/bin/env python3
"""This script filters configurations with large forces."""
from ...mlip.get_properties import get_max_force
from a_parameter import (input_config, atom_symbols_in_cfg,
                         force_length_criterion)
# , removed_cfg_file,
#                          kept_cfg_file)

if __name__ == '__main__':

    print(f'File: {input_config}')
    get_max_force(input_config, atom_symbols_in_cfg, 'fig-force.pdf',
                  length_or_component='length',
                  force_criterion=force_length_criterion)
                  # removed_cfg_file=removed_cfg_file,
                  # kept_cfg_file=kept_cfg_file)
