"""This script contains functions for reading SiO2 configurations."""
import numpy as np
from ase.io import read
from ..util.util import set_actual_atom_symbols
from ..util.SiO2_parameter import (Si_O_H_Al_atom_symbol_tuple_lammps, atomic_number_to_POTCAR_order)


def read_SiO2_dump(file, index):
    """This function reads a lammps dump file of SiO2 structure and changes the
    atom index correctly. A list of Atoms-type objects is returned."""
    # 'configs' is a list of Atoms-type objects.
    configs = read(file, format='lammps-dump-text', index=index)
    if type(configs).__name__ == 'Atoms':
        configs = [configs]
        # If there is one configuration, convert to a list of configs.
    for config in configs:
        config = set_actual_atom_symbols(config,
                Si_O_H_Al_atom_symbol_tuple_lammps)
    return configs


def sort_config_by_POTCAR_order(config, return_sort_index=False):
    """This function sorts atoms in a configuration by the POTCAR order.
    When 'return_sort_index' == True, sort_index is also returned."""
    number_to_sort = [atomic_number_to_POTCAR_order[atomic_number]
            for atomic_number in config.get_atomic_numbers()]
    # number_to_sort = np.zeroes(len(atomic_numbers))
    # for index, atomic_number in enumerate(atomic_numbers):
    #     number_to_sort[index] = atomic_number_to_POTCAR_order[atomic_number]
    sort_index = np.argsort(number_to_sort)
    config = config[sort_index]
    if return_sort_index:
        return_value = config, sort_index
    else:
        return_value = config
    return return_value
