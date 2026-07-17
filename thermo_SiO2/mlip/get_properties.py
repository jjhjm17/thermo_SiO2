"""This script obtains properties from configurations from SiO2 configurations."""
from sys import argv
import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
from ase.io import read
from .read import read_SiO2
from ..mlip.write_config import write_cfg_SiO2


def get_max_force(file, atom_symbols, fig_file='fig-force-hist.pdf', *,
                  length_or_component, force_criterion=None,
                  removed_cfg_file=None, kept_cfg_file=None,
                  std_or_max):
    """This function obtains the maximum force. This information can be used for the filtering.
    std_or_max: 'std' or 'max', watch the standard deviation of forces or
        the maximum force."""
    VERBOSE = False
    if std_or_max == 'std' and length_or_component == 'length':
        print("Error: for std force, only length_or_component='component' is "
              "supported, because we think mean is near 0.")
    if 'OUTCAR' in file or 'vasprun.xml' in file:
        configs = read(file, index=':')
    else:
        configs = read_SiO2(file, atom_symbols)
    # global_max_force = 0
    global_watch_force = 0  # watch_force is std or max force

    # all_forces = np.array([])
    # max_forces = np.array([])  # max forces for each configuration.
    watch_forces = np.array([])  # max forces for each configuration.
    kept_forces = np.array([])  # max forces for each configuration.
    print()
    kept_cfgs = []
    removed_cfgs = []
    print(f'std_or_max force = {std_or_max}')
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
        # all_forces = np.concatenate((all_forces, forces))
        if std_or_max == 'std':
            watch_force_in_config = np.std(forces)
            # mean will be near zero
        elif std_or_max == 'max':
            # max_force_in_config = np.max(forces)
            forces = np.absolute(forces)
            watch_force_in_config = np.max(forces)
        if force_criterion is not None:
            if watch_force_in_config > force_criterion:
                if VERBOSE:
                    print(f'\n Configuration with the {std_or_max}. force larger than the '
                        'criterion:')
                    print(f'Index = {index}')
                    print(f'{config = }')
                    print(f'watched (max or std) force = {watch_force_in_config:.6g} eV/Ang\n')
                    print(f'{config.get_cell() = }')
                    print(f'{config.get_positions()[:5] = }')
                removed_cfgs.append(config)
            else:
                kept_cfgs.append(config)
                kept_forces = np.append(kept_forces, watch_force_in_config)
        watch_forces = np.append(watch_forces, watch_force_in_config)

    if removed_cfgs != None:
        print(f'\n number of initial cfgs: {len(configs)}')
        print(f'removed cfgs: {len(removed_cfgs)}')
        print(f'kept cfgs: {len(kept_cfgs)}\n')
        write_cfg_SiO2(removed_cfg_file, removed_cfgs)
        write_cfg_SiO2(kept_cfg_file, kept_cfgs)
    
        print(f'Removed cfgs saved: {removed_cfg_file}')
        print(f'Kept cfgs saved: {kept_cfg_file}')

    global_watch_force = np.max(watch_forces)
    print(f'\n Global max {std_or_max} force = {global_watch_force:.6g} eV/Ang')

    # all_forces = np.array(all_forces).flatten()
    fig, axes = plt.subplots(2, 1, tight_layout=True)
    # ax.hist(all_forces, bins=30)
    ax = axes[0]
    ax.hist(watch_forces, bins=30)
    ax.set_ylabel('Frequency of initial cfgs')
    # ax.set_xlabel('Force components (eV/Ang)')

    ax = axes[1]
    ax.hist(kept_forces, bins=30)
    ax.set_ylabel('Frequency of kept cfgs')
    ax.set_xlabel(f'{std_or_max}. force in a snapshot (eV/Ang)')
    plt.savefig(fig_file, transparent=True)
    print(f'File saved: {fig_file}')
    plt.show()

    return global_watch_force
