#!/usr/bin/env python3
from ...util.util import shell

# Ref
# https://github.com/skelton-group/Phonopy-Spectroscopy/tree/master/example/a-SiO2

shell('phonopy --dim="6 6 3" --readfc --hdf5 --fc-symmetry --mesh="1 1 1" --eigenvectors')

