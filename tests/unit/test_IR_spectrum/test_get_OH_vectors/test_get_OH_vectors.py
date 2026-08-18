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
import os
import tempfile
import unittest
import numpy as np
import yaml
from ase import Atoms
from thermo_SiO2.io import read_sil
from thermo_SiO2.IR_spectrum.dipole.get_OH_vectors import (
    classify_OH_types,
    get_OH_vectors,
    find_bonded_oxygens,
    select_H_indices,
)


class TestGetOHVectors(unittest.TestCase):

    @staticmethod
    def load_OH_analysis(path):
        rows = []
        with open(path) as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                h_idx, OH_type, *numeric_values = line.split()
                rows.append(
                    (int(h_idx), OH_type)
                    + tuple(float(value) for value in numeric_values)
                )
        return rows

    def test_classify_OH_types(self):
        cases = {
            'Si-OH': ['Si'],
            'Al-OH': ['Al'],
            'Al-(OH)-Si': ['Al', 'Si'],
            'Al-(OH)-Al': ['Al', 'Al'],
            'Si-(OH)-Si': ['Si', 'Si'],
            'other': [],
        }

        for expected, cations in cases.items():
            with self.subTest(expected=expected):
                symbols = ['O', 'H'] + cations
                positions = [[5.0, 5.0, 5.0], [5.0, 5.9, 5.0]]
                cation_positions = [
                    [3.4, 5.0, 5.0],
                    [6.6, 5.0, 5.0],
                ]
                positions.extend(cation_positions[:len(cations)])
                atoms = Atoms(
                    symbols=symbols,
                    positions=positions,
                    cell=[20.0, 20.0, 20.0],
                    pbc=True,
                )
                self.assertEqual(
                    classify_OH_types(atoms, [0]), [expected]
                )

    def test_select_H_indices(self):
        atoms = Atoms('OHHSiH')

        self.assertEqual(
            select_H_indices({'OH_vector_H_indices': [2]}, atoms),
            [2],
        )
        self.assertEqual(
            select_H_indices({'OH_vector_H_all': True}, atoms),
            [1, 2, 4],
        )

        with self.assertRaisesRegex(ValueError, 'not both'):
            select_H_indices(
                {
                    'OH_vector_H_indices': [2],
                    'OH_vector_H_all': True,
                },
                atoms,
            )
        with self.assertRaisesRegex(ValueError, 'Set OH_vector_H_indices'):
            select_H_indices({}, atoms)

    def test_classify_H2O_requires_isolated_H_O_H_bonds(self):
        water = Atoms(
            symbols=['O', 'H', 'H'],
            positions=[
                [5.0, 5.0, 5.0],
                [5.0, 5.9, 5.0],
                [5.0, 4.1, 5.0],
            ],
            cell=[20.0, 20.0, 20.0],
            pbc=True,
        )
        self.assertEqual(classify_OH_types(water, [0]), ['H2O'])

        water_with_extra_H_bond = water.copy()
        water_with_extra_H_bond += Atoms(
            symbols=['Si'], positions=[[5.0, 7.0, 5.0]]
        )
        self.assertNotEqual(
            classify_OH_types(water_with_extra_H_bond, [0]), ['H2O']
        )

    def test_classify_Al_OH2(self):
        atoms = Atoms(
            symbols=['O', 'H', 'H', 'Al'],
            positions=[
                [5.0, 5.0, 5.0],
                [5.0, 5.9, 5.0],
                [5.0, 4.0, 5.0],
                [7.0, 5.0, 5.0],
            ],
            cell=[20.0, 20.0, 20.0],
            pbc=True,
        )
        self.assertEqual(classify_OH_types(atoms, [0]), ['Al..OH2'])

    def test_classify_OH_types_across_periodic_boundary(self):
        atoms = Atoms(
            symbols=['O', 'H', 'Si'],
            positions=[
                [0.2, 5.0, 5.0],
                [0.2, 5.9, 5.0],
                [9.0, 5.0, 5.0],
            ],
            cell=[10.0, 10.0, 10.0],
            pbc=True,
        )
        self.assertEqual(classify_OH_types(atoms, [0]), ['Si-OH'])

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

    def test_write_OH_analysis(self):
        with open('fixtures/in.yaml') as f:
            param = yaml.safe_load(f)

        with tempfile.TemporaryDirectory() as tmp_dir:
            analysis_file = os.path.join(tmp_dir, 'OH_analysis.dat')
            in_file = os.path.join(tmp_dir, 'in.yaml')
            param['OH_analysis_out'] = analysis_file
            with open(in_file, 'w') as f:
                yaml.safe_dump(param, f)

            get_OH_vectors(in_file)

            with open(analysis_file) as f:
                lines = f.read().splitlines()

        self.assertEqual(
            lines[0],
            '# H_index OH_type OH_bond_length (Ang) '
            'Wrapped_x Wrapped_y Wrapped_z (Ang)',
        )
        self.assertEqual(lines[1].split()[:2], ['1', 'other'])
        self.assertEqual(lines[2].split()[:2], ['3', 'other'])
        self.assertAlmostEqual(float(lines[1].split()[2]), 1.0)
        self.assertAlmostEqual(float(lines[2].split()[2]), 1.0)
        np.testing.assert_allclose(
            [float(value) for value in lines[1].split()[3:]],
            [2.0, 1.0, 1.0],
            atol=1e-5,
        )

    def test_OH_analysis_against_hand_fixture(self):
        fixture_dir = 'fixtures_b_OH_analysis'
        with open(os.path.join(fixture_dir, 'in.yaml')) as f:
            param = yaml.safe_load(f)

        with tempfile.TemporaryDirectory() as tmp_dir:
            analysis_file = os.path.join(tmp_dir, 'OH_analysis.dat')
            in_file = os.path.join(tmp_dir, 'in.yaml')
            param['OH_analysis_out'] = analysis_file
            with open(in_file, 'w') as f:
                yaml.safe_dump(param, f)

            get_OH_vectors(in_file)
            actual = self.load_OH_analysis(analysis_file)

        expected = self.load_OH_analysis(
            os.path.join(fixture_dir, 'OH_analysis_hand.dat')
        )
        self.assertEqual(
            [(row[0], row[1]) for row in actual],
            [(row[0], row[1]) for row in expected],
        )
        np.testing.assert_allclose(
            [row[2:] for row in actual],
            [row[2:] for row in expected],
            atol=1e-5,
            rtol=0,
        )


if __name__ == '__main__':
    unittest.main()
