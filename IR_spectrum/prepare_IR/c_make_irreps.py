#!/usr/bin/env python3
from parameters import supercell_dim
from ...util.util import shell

# Ref
# https://github.com/skelton-group/Phonopy-Spectroscopy/tree/master/example/a-SiO2

shell(f'phonopy --dim="{supercell_dim}" --readfc --hdf5 --fc-symmetry --irreps="0 0 0"')

