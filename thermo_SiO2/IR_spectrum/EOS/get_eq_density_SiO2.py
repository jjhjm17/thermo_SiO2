#!/usr/bin/env python3
"""This script obtains the density of the equilibrium structure from the energy volume curve."""
import numpy as np


def get_eq_density_SiO2(EVinet_path, approx_eq_config):
    """This function obtains equilibrium density."""
    eq_vol = np.loadtxt(EVinet_path)[1]  # Ang3/atom
    # eq_vol = 14.492
    # print(f'equilibrium volume of Vinet curve = {eq_vol:.5g} Ang^3/atom')
    formula = approx_eq_config.get_chemical_formula(empirical=True)
    if formula != 'O2Si':
        print('Error: please implement for non SiO2 structure.')
    mass_Si = 28.0855  # getAtomicMass.sh Si
    mass_O = 15.9994  # getAtomicMass.sh O
    # N_A = 6.02214076e23  # Avogadro's constant
    # density = (mass_Si + 2 * mass_O) / 3 / (eq_vol * (1e-9)**3 * N_A) *   # g/m^3
    density = (mass_Si + 2 * mass_O) / 3 / eq_vol * 1.6605391  # g/(cm)^3
    # units program
    # You have: amu/(angstrom^3)
    # You want: g/(cm)^3
    print(f'density = {density:.5g} g/cm^3')
