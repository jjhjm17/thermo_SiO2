#!/usr/bin/env python3
import os
from ase.io import read
from parameters import supercell_dim
from ...util.util import shell

# anime-tag
# Gamma-point animation
# shell(f'phonopy --hdf5 --readfc --dim="{supercell_dim}" --anime="0 0 0" ')

# Ref. https://henriquemiranda.github.io/phononwebsite/index.html
shell(f'phonopy --fc-format=hdf5 --readfc --dim="{supercell_dim}"  --band-format=yaml  band.conf')
