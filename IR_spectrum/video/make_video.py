#!/usr/bin/env python3
import os
import shutil
from ase.io import read
from parameters import supercell_dim, sample_folder
from ...util.util import shell
import parameters

# anime-tag
# Gamma-point animation
# shell(f'phonopy --hdf5 --readfc --dim="{supercell_dim}" --anime="0 0 0" ')

shutil.copy(f'{sample_folder}/POSCAR', '.')
os.symlink(f'{sample_folder}/force_constants.hdf5', 'force_constants.hdf5')

if hasattr(parameters, 'format'):
    format = parameters.format
else:
    format = 'TSS_physics_website'

if format == 'jmol':
    shell(f'phonopy --fc-format=hdf5 --readfc --dim="{supercell_dim}"  setting_jmol.conf')

elif format == 'TSS_physics_website':
    # Ref. https://henriquemiranda.github.io/phononwebsite/index.html
    shell(f'phonopy --fc-format=hdf5 --readfc --dim="{supercell_dim}"  --band-format=yaml  band.conf')
