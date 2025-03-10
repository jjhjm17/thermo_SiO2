"""This script obtains properties from configurations from SiO2 configurations."""
from sys import argv
import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
from ase.io import read
from .read import read_SiO2
from ..mlip.write_config import write_cfg_SiO2


def get_max_force(file, atom_symbols, fig_file='fig-force-hist.pdf',
                  length_or_component='length', force_criterion=None):
                  # removed_cfg_file=None, kept_cfg_file=None):
    """This function obtains the maximum force. This information can be used for the filtering."""
    if 'OUTCAR' in file or 'vasprun.xml' in file:
        configs = read(file, index=':')
    else:
        configs = read_SiO2(file, atom_symbols)
    global_max_force = 0

    # all_forces = np.array([])
    max_forces = np.array([])  # max forces for each configuration.
    print()
    kept_cfgs = []
    removed_cfgs = []
    for index, config in enumerate(configs):
        if hasattr(config, 'forces'):
            forces = config.forces
        else:
            forces = config.get_forces()
        if length_or_component == 'length':
            force_lengths = [norm(vector) for vector in forces]
            forces = force_lengths
        else:
            forces = forces.flatten()
        forces = np.absolute(forces)
        # all_forces = np.concatenate((all_forces, forces))
        max_force_in_config = np.max(forces)
        if force_criterion is not None:
            if max_force_in_config > force_criterion:
                print('\n Configuration with the max. force larger than the '
                      'criterion:')
                print(f'Index = {index}')
                print(f'{config = }')
                print(f'max. force = {max_force_in_config:.6g} eV/Ang\n')
                print(f'{config.get_cell() = }')
                print(f'{config.get_positions()[:5] = }')
                removed_cfgs.append(config)
            else:
                kept_cfgs.append(config)
        max_forces = np.append(max_forces, max_force_in_config)

    # write_cfg_SiO2(removed_cfg_file, removed_cfgs)
    # write_cfg_SiO2(kept_cfg_file, kept_cfgs)
    
    # print(f'Removed cfgs saved: {removed_cfg_file}')
    # print(f'Kept cfgs saved: {kept_cfg_file}')

    global_max_force = np.max(max_forces)
    print(f'\n Global max force = {global_max_force:.6g} eV/Ang')

    # all_forces = np.array(all_forces).flatten()
    fig, ax = plt.subplots(1, 1, tight_layout=True)
    # ax.hist(all_forces, bins=30)
    ax.hist(max_forces, bins=30)
    # ax.set_xlabel('Force components (eV/Ang)')
    ax.set_xlabel('Norm of max. force in a snapshot (eV/Ang)')
    ax.set_ylabel('Frequency')
    plt.savefig(fig_file, transparent=True)
    print(f'File saved: {fig_file}')
    plt.show()

    return global_max_force
