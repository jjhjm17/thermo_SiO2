#!/usr/bin/env python3
import os
from parameters import supercell_dim
from ...util.util import shell

# Ref
# https://github.com/skelton-group/Phonopy-Spectroscopy/tree/master/example/a-SiO2

# command = 'phonopy --dim=' f'"{supercell_dim}"' ' --readfc --hdf5 --fc-symmetry --mesh="1 1 1" --eigenvectors'
command = 'phonopy  --readfc-format hdf5 --hdf5 --no-sym-fc --mesh="1 1 1" --eigenvectors'
# print(f'{command=}')
shell(command)

