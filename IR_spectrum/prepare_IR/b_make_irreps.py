#!/usr/bin/env python3
from ...util.util import shell

# Ref
# https://github.com/skelton-group/Phonopy-Spectroscopy/tree/master/example/a-SiO2

shell('phonopy --dim="6 6 3" --readfc --hdf5 --fc-symmetry --irreps="0 0 0"')

