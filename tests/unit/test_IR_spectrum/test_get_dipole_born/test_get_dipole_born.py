import unittest
import numpy as np
from thermo_SiO2.IR_spectrum.dipole.get_dipole_born import get_dipole_born

class TestGetDipoleBorn(unittest.TestCase):

    def test_get_dipole_born_fixed_q(self):

        dipoles = get_dipole_born('fixtures_a_fixed_q/in.yaml')
        dipole_lammps = np.loadtxt('fixtures_a_fixed_q/dipole_lammps.dat')
        diff = dipole_lammps[:, 1:] - dipoles 
        # align 1st values
        diff[:] -= dipole_lammps[0, 1:] - dipoles[0]
        diff_norm = np.linalg.norm(diff) / len(diff)
        self.assertLess(diff_norm, 1e-3, "Should be small")

    def test_get_dipole_born(self):

        dipoles = get_dipole_born('fixtures_b_born/in.yaml')
        dipole_hand = np.loadtxt('fixtures_b_born/dipole_hand.dat')
        # hand + claude
        diff = dipole_hand[:, 1:] - dipoles 
        diff_norm = np.linalg.norm(diff) / len(diff)
        self.assertLess(diff_norm, 1e-3, "Should be small")

if __name__ == '__main__':
    unittest.main()

