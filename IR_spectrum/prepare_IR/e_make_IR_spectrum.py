#!/usr/bin/env python3
from ...util.util import shell
from parameters import temperature, freq_range, linewidth

# Ref
# https://github.com/skelton-group/Phonopy-Spectroscopy/tree/master/example/a-SiO2

# shell(f'phonopy-ir --linewidth-hdf5="{linewidth_file}" --linewidth-temperature={temperature} --spectrum-range="{freq_range}"')  # range: cm^-1
shell(f'phonopy-ir --linewidth="{linewidth}" --spectrum-range="{freq_range}"')  # range: cm^-1
