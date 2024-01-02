#!/usr/bin/env python3
from ...util.util import shell
from parameters import freq_range

# Ref
# https://github.com/skelton-group/Phonopy-Spectroscopy/tree/master/example/a-SiO2

shell(f'phonopy-ir --linewidth-hdf5="kappa-m484848-g0.hdf5" --linewidth-temperature=300 --spectrum-range="{freq_range}"')  # range: cm^-1
