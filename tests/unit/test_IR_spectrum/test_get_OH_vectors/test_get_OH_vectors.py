"""
This script tests get_OH_vectors for a simple case.
Input files: In ./fixtures/ are there files.
  in.yaml: input file to invoke getting OH vectors
  config.dump:
    10 x 10 x 5 Ang, hexagonal cell.
    4 atoms (two OH's).
       At step 0
         O1 is at (1,1,1), H1 at (2,1,1)
         O2 is at (3,3,4), H2 at (3,3,10)  # shift is nonzero
       At step 1
         O1 is at (1.1,1.1,1.1), H1 at (2.2,1.2,1.2)
         O2 is at (3.3,3.3,4.3), H2 at (3.5,3.5,10.5)
  OH_vectors_hand.dat : It has OH_vectors to be compared.
Draft written by Claude, modified by hands.
"""
import unittest
import numpy as np
from thermo_SiO2.io import read_sil
from thermo_SiO2.IR_spectrum.dipole.get_OH_vectors import (
    get_OH_vectors,
    find_bonded_oxygens,
)


class TestGetOHVectors(unittest.TestCase):

    def test_get_OH_vectors(self):

        print('\ntest_get_OH_vectors')
        # oh_vectors = np.asarray(get_OH_vectors('fixtures/in.yaml'))
        oh_vectors = get_OH_vectors('fixtures/in.yaml')
        # oh_vectors shape: (n_frames, n_pairs, 3) -> flatten to (n_frames, n_pairs*3)
        n_frames, n_pairs, _ = oh_vectors.shape
        oh_vectors_flat = oh_vectors.reshape(n_frames, n_pairs * 3)

        oh_hand = np.loadtxt('fixtures/OH_vectors_hand.dat')
        diff = oh_hand[:, 1:] - oh_vectors_flat
        diff_norm = np.linalg.norm(diff) / len(diff)
        self.assertLess(diff_norm, 1e-3, "Should be small")

    def test_find_bonded_oxygens(self):
        # H index 1 (H1) should bond to O index 0 (O1);
        # H index 3 (H2) should bond to O index 2 (O2).
        # Both bonds are exactly 1.0 Ang, well under the 1.5 Ang cutoff,
        # and each H has only one O within range.
        print('\ntest_find_bonded_oxygens')
        cfgs = read_sil('fixtures/config.dump', atom_symbols='Si O H Al')
        cfg0 = cfgs[0]

        H_indices = [1, 3]
        O_indices, shifts = find_bonded_oxygens(cfg0, H_indices)

        self.assertEqual(O_indices, [0, 2])
        self.assertEqual(shifts.shape, (2, 3))

        # The shift of 2nd H is nonzero.
        np.testing.assert_allclose(shifts, np.array([[0, 0, 0], [0, 0, -5]]),
                                   atol=1e-8)

        # # Double check the bond lengths found are close to 1.0 Ang.
        # positions0 = cfg0.get_positions()
        # for h_idx, o_idx in zip(H_indices, O_indices):
        #     bond_len = np.linalg.norm(positions0[h_idx] - positions0[o_idx])
        #     breakpoint()
        #     # np.testing.assert_allclose(bond_len, 1, 0.8)


if __name__ == '__main__':
    unittest.main()
