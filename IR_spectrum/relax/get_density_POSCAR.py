#!/usr/bin/env python3
"""This script postprocesses: obtains densities of SiO2, get POSCAR, and extrapolation grade."""

import os
import numpy as np
from ase.io import read, write
from mlippy.atms import ase_savecfgs
from a_parameters import calc_folder, structures
from ...util.util import read_output, set_actual_atom_symbols
from ...util.SiO2_parameter import Si_O_H_Al_atom_symbol_tuple_lammps
from ...mlip.read_config import sort_config_by_POTCAR_order
from ...mlip.write_config import write_cfg_SiO2


def get_density_POSCAR():
    """This function cleans calculation files."""
    print(f'{structures = }\n')

    os.chdir(calc_folder)
    densities = []
    print('Folder   number of configs')
    # for index_seed in range(1, num_seeds+1):  # 1 to num_seeds
    for i_config, config_file in enumerate(structures):
        calc_subfolder = f'struct_{str(i_config).zfill(3)}'
        os.chdir(calc_subfolder)
        # num_configs = read_output('grep TIMESTEP dump.atom | wc --line')
        # print(f'{calc_subfolder} {num_configs}')
        print(f'{calc_subfolder} 1')
        config = read('dump.atom', index='-1', format='lammps-dump-text')
        config = set_actual_atom_symbols(config,
            Si_O_H_Al_atom_symbol_tuple_lammps)
        amu_per_ang_3__to__g_per_cm_3 = 1.6605391
        # units u/angstrom^3 g/cm^3
        density = (np.sum(config.get_masses()) / config.get_volume() *
                amu_per_ang_3__to__g_per_cm_3)
        densities.append(density)
        config = sort_config_by_POTCAR_order(config)
        print(f'{config=}')
        write('POSCAR', config, format='vasp', vasp5=True)
        write_cfg_SiO2('relaxed.cfg', [config],
                       desc='SiOAl melt quench by Julian PBE, relaxed by vdW MTP')

        os.chdir('..')
    os.chdir('..')
    avg_density = np.average(densities)
    std = np.std(densities, ddof=1)
    print(f'Density = {avg_density:.4f} +- {std:.4f} (std) g/cm^3 '
          f'(n={len(densities)})')
    print("written: POSCAR, relaxed.cfg.") 

if __name__ == "__main__":
    get_density_POSCAR()
