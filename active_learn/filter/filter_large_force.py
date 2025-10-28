#!/usr/bin/env python3
"""This script filters configurations with large forces."""
from ...mlip.get_properties import get_max_force
from a_parameter import (input_config, atom_symbols_in_cfg,
                         force_length_criterion, std_or_max, length_or_component)
import a_parameter as param
# , removed_cfg_file,
#                          kept_cfg_file)

if __name__ == '__main__':

    print(f'File: {input_config}')
    if hasattr(param, 'removed_cfg_file'):
        removed_cfg_file = param.removed_cfg_file
    else:
        removed_cfg_file = None
    if hasattr(param, 'kept_cfg_file'):
        kept_cfg_file = param.kept_cfg_file
    else:
        kept_cfg_file = None
    get_max_force(input_config, atom_symbols_in_cfg, 'fig-force.pdf',
                  length_or_component=length_or_component,
                  force_criterion=force_length_criterion,
                  removed_cfg_file=removed_cfg_file,
                  kept_cfg_file=kept_cfg_file,
                  std_or_max=std_or_max)
