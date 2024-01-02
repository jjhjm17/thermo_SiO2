#!/usr/bin/env python3
from ...util.util import shell

# Ref
# https://github.com/skelton-group/Phonopy-Spectroscopy/tree/master/example/a-SiO2

shell('phono3py --dim="2 2 2" --dim-fc2="6 6 3" --fc2 --fc3 -v --br --thm --mesh="48 48 48" --write-gamma --gp=0')
