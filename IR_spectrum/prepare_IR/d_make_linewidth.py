#!/usr/bin/env python3
from parameters import supercell_dim, supercell_dim_fc3, mesh, temperature
from ...util.util import shell

# Ref
# https://github.com/skelton-group/Phonopy-Spectroscopy/tree/master/example/a-SiO2

shell(
    # f'phono3py --dim="{supercell_dim_fc3}" --dim-fc2="{supercell_dim}" --fc2 --fc3 -v --br --thm --mesh="{mesh}" --write-gamma --gp=0 --ts="{temperature}"'
    f'phono3py --dim="{supercell_dim_fc3}" --dim-fc2="{supercell_dim}" --fc2 --fc3 -v --br --sigma="0.1" --mesh="{mesh}" --write-gamma --gp=0 --ts="{temperature}"'
)

# --br (BTERTA = .TRUE.)
# Run calculation of lattice thermal conductivity tensor with the single mode relaxation time approximation (RTA) and linearized phonon Boltzmann equation.
# A direct solution requires large memory.

# --write-gamma (WRITE_GAMMA = .TRUE.)
# Imaginary parts of self energy at harmonic phonon frequencies
#  are written into file in hdf5 format. The result is written into kappa-mxxx-gx(-sx-sdx).hdf5 or kappa-mxxx-gx-bx(-sx-sdx).hdf5

# --gp (GRID_POINTS)
# Grid points are specified by their unique indices, e.g., for selecting the q-points where imaginary parts of self energees are calculated.

# https://phonopy.github.io/phono3py/command-options.html#br-bterta-true
