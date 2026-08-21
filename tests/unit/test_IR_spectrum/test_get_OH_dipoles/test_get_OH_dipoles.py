"""
This script tests get_OH_dipoles for a simple case.
Input files: In ./fixtures/ are there files.
  in.yaml: input file to invoke getting OH dipoles
  config.dump:
    10 x 10 x 5 Ang, hexagonal cell.
    4 atoms (two OH's).
       At step 0
         O1 is at (1,1,1), H1 at (2,1,1)
         O2 is at (3,3,4), H2 at (3,3,10)  # shift is nonzero
       At step 1
         O1 is at (1.1,1.1,1.1), H1 at (2.2,1.2,1.2)
         O2 is at (3.3,3.3,4.3), H2 at (3.5,3.5,10.5)
  OH_dipoles_hand.dat : It has OH_dipoles to be compared.
Draft written by Claude, modified by hands.
"""
import json
import os
import tempfile
import unittest
import numpy as np
import pandas as pd
import yaml
from ase import Atoms
from ase.io import write
from thermo_SiO2.io import read_sil
from thermo_SiO2.IR_spectrum.dipole.get_OH_dipoles import (
    build_dipole_groups,
    classify_OH_types,
    compute_frame_dipoles,
    find_hydrogen_bonds,
    get_OH_dipoles,
    get_charge_tensors,
    find_bonded_oxygens,
    select_H_indices,
    write_OH_analysis,
)


class TestGetOHDipoles(unittest.TestCase):

    @staticmethod
    def load_OH_analysis(path):
        analysis = pd.read_csv(path)
        rows = analysis.to_dict(orient='records')
        for row in rows:
            row['Hbond_lengths'] = json.loads(row['Hbond_lengths'])
            row['Hbond_O_indices'] = json.loads(row['Hbond_O_indices'])
        return rows

    def assert_OH_dipoles_match_fixture(self, fixture_dir):
        """Compare calculated OH dipoles with a hand-computed fixture."""
        oh_dipoles = get_OH_dipoles(os.path.join(fixture_dir, 'in.yaml'))
        oh_hand = np.loadtxt(
            os.path.join(fixture_dir, 'OH_dipoles_hand.dat')
        )

        n_frames, n_groups, n_components = oh_dipoles.shape
        self.assertEqual(n_components, 3)
        self.assertEqual(
            oh_hand.shape,
            (n_frames, 1 + n_groups * n_components),
        )
        np.testing.assert_array_equal(oh_hand[:, 0], np.arange(n_frames))
        np.testing.assert_allclose(
            oh_hand[:, 1:],
            oh_dipoles.reshape(n_frames, n_groups * n_components),
            atol=1e-5,
            rtol=0,
        )

    def test_find_hydrogen_bonds(self):
        atoms = Atoms(
            symbols=['O', 'H', 'O', 'O', 'O', 'O'],
            positions=[
                [5.0, 5.0, 5.0],
                [6.0, 5.0, 5.0],
                [7.5, 5.0, 5.0],
                [8.2, 5.0, 5.0],
                [6.0, 7.0, 5.0],
                [8.6, 5.0, 5.0],
            ],
            cell=[20.0, 20.0, 20.0],
            pbc=True,
        )
        matches = find_hydrogen_bonds(atoms, [1], [0])[0]
        self.assertEqual([O_idx for O_idx, _ in matches], [2, 3])
        np.testing.assert_allclose(
            [distance for _, distance in matches], [1.5, 2.2]
        )
        # O index 5 has H..O = 2.6 Ang, but donor-O..acceptor-O =
        # 3.6 Ang, so the corrected O..O criterion excludes it.

        no_matches = find_hydrogen_bonds(atoms, [1], [0], min_angle_deg=181)
        self.assertEqual(no_matches, [[]])

        with tempfile.TemporaryDirectory() as tmp_dir:
            analysis_file = os.path.join(tmp_dir, 'OH_analysis.csv')
            write_OH_analysis(
                analysis_file,
                [1],
                ['Si-OH'],
                [1.0],
                [matches],
                [[6.0, 5.0, 5.0]],
            )
            analysis = pd.read_csv(analysis_file)
        self.assertEqual(
            json.loads(analysis.loc[0, 'Hbond_lengths']), [1.5, 2.2]
        )
        self.assertEqual(
            json.loads(analysis.loc[0, 'Hbond_O_indices']), [2, 3]
        )

    def test_find_hydrogen_bond_across_periodic_boundary(self):
        atoms = Atoms(
            symbols=['O', 'H', 'O'],
            positions=[
                [8.5, 5.0, 5.0],
                [9.5, 5.0, 5.0],
                [1.0, 5.0, 5.0],
            ],
            cell=[10.0, 10.0, 10.0],
            pbc=True,
        )
        matches = find_hydrogen_bonds(atoms, [1], [0])[0]
        self.assertEqual(matches[0][0], 2)
        self.assertAlmostEqual(matches[0][1], 1.5)

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
            select_H_indices({'OH_dipole_H_indices': [2]}, atoms),
            [2],
        )
        self.assertEqual(
            select_H_indices({'OH_dipole_H_all': True}, atoms),
            [1, 2, 4],
        )

        with self.assertRaisesRegex(ValueError, 'not both'):
            select_H_indices(
                {
                    'OH_dipole_H_indices': [2],
                    'OH_dipole_H_all': True,
                },
                atoms,
            )
        with self.assertRaisesRegex(ValueError, 'Set OH_dipole_H_indices'):
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

    def test_two_H_groups_use_smaller_index_and_sum_dipoles(self):
        cases = [
            (
                'H2O',
                ['O', 'H', 'H'],
                [[5.0, 5.0, 5.0], [5.8, 5.5, 5.0], [4.6, 4.2, 5.0]],
            ),
            (
                'Al..OH2',
                ['O', 'H', 'H', 'Al'],
                [
                    [5.0, 5.0, 5.0],
                    [5.8, 5.5, 5.0],
                    [4.6, 4.2, 5.0],
                    [7.0, 5.0, 5.0],
                ],
            ),
        ]
        for expected_type, symbols, positions in cases:
            with self.subTest(expected_type=expected_type):
                atoms = Atoms(
                    symbols=symbols,
                    positions=positions,
                    cell=[20.0, 20.0, 20.0],
                    pbc=True,
                )
                groups, pair_by_H = build_dipole_groups(atoms, [2, 1])
                self.assertEqual(len(groups), 1)
                self.assertEqual(groups[0]['H_index'], 1)
                self.assertEqual(groups[0]['H_indices'], [1, 2])
                self.assertEqual(groups[0]['OH_type'], expected_type)

                charge_tensors = np.repeat(
                    np.eye(3)[None, :, :], len(atoms), axis=0
                )
                actual = compute_frame_dipoles(
                    atoms.positions, groups, pair_by_H, charge_tensors
                )[0]
                expected = (
                    atoms.positions[1] - atoms.positions[0]
                    + atoms.positions[2] - atoms.positions[0]
                )
                np.testing.assert_allclose(actual, expected)

    def test_charge_tensors(self):
        atoms = Atoms(
            symbols=['O', 'H', 'H'],
            positions=[[5.0, 5.0, 5.0], [5.8, 5.5, 5.0], [4.6, 4.2, 5.0]],
            cell=[20.0, 20.0, 20.0],
            pbc=True,
        )
        H_tensor = np.array([
            [2.0, 1.0, 0.0],
            [0.0, 4.0, 0.0],
            [0.0, 0.0, 6.0],
        ])
        with tempfile.TemporaryDirectory() as tmp_dir:
            poscar = os.path.join(tmp_dir, 'POSCAR')
            born = os.path.join(tmp_dir, 'BORN')
            write(poscar, atoms, format='vasp', direct=True, vasp5=True)
            with open(born, 'w') as f:
                f.write(
                    '# epsilon and Z*\n'
                    '1 0 0 0 1 0 0 0 1\n'
                    '-2 0 0 0 -2 0 0 0 -2\n'
                    '2 1 0 0 4 0 0 0 6\n'
                    '3 0 0 0 3 0 0 0 3\n'
                )

            nominal = get_charge_tensors({'charge': 'nominal'}, atoms)
            isotropic = get_charge_tensors(
                {
                    'charge': 'born_isotropic',
                    'born_file': born,
                    'born_poscar': poscar,
                },
                atoms,
            )
            full = get_charge_tensors(
                {
                    'charge': 'born_full',
                    'born_file': born,
                    'born_poscar': poscar,
                },
                atoms,
            )

        np.testing.assert_allclose(nominal[1], np.eye(3))
        np.testing.assert_allclose(isotropic[1], 4.0 * np.eye(3))
        np.testing.assert_allclose(full[1], H_tensor)
        with self.assertRaisesRegex(ValueError, "Missing required 'charge'"):
            get_charge_tensors({}, atoms)
        with self.assertRaisesRegex(ValueError, 'Invalid charge'):
            get_charge_tensors({'charge': 'unknown'}, atoms)
        with self.assertRaisesRegex(ValueError, 'requires born_file'):
            get_charge_tensors({'charge': 'born_full'}, atoms)


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

    def test_get_OH_dipoles(self):
        print('\ntest_get_OH_dipoles')
        self.assert_OH_dipoles_match_fixture('fixtures')

    def test_get_OH_dipoles_born_isotropic(self):
        self.assert_OH_dipoles_match_fixture(
            'fixtures_c_born_isotropic'
        )

    def test_get_OH_dipoles_born_full(self):
        self.assert_OH_dipoles_match_fixture('fixtures_d_born_full')

    def test_write_BORN_isotropic(self):
        fixture_dir = 'fixtures_c_born_isotropic'
        with open(os.path.join(fixture_dir, 'in.yaml')) as f:
            param = yaml.safe_load(f)

        with tempfile.TemporaryDirectory() as tmp_dir:
            born_isotropic_file = os.path.join(
                tmp_dir, 'BORN_isotropic.txt'
            )
            in_file = os.path.join(tmp_dir, 'in.yaml')
            param['BORN_isotropic_out'] = born_isotropic_file
            with open(in_file, 'w') as f:
                yaml.safe_dump(param, f)

            get_OH_dipoles(in_file)

            with open(born_isotropic_file) as f:
                header = f.readline().rstrip()
            actual = np.loadtxt(born_isotropic_file)

        self.assertEqual(header, '# H_index Born isotropic charge')
        np.testing.assert_allclose(
            actual,
            [[1, 1.0], [3, 2.0]],
            atol=1e-5,
            rtol=0,
        )

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

    def test_find_bonded_oxygens_selects_nearest_with_warning(self):
        atoms = Atoms(
            symbols=['O', 'O', 'H'],
            positions=[
                [0.0, 0.0, 0.0],
                [2.48, 0.0, 0.0],
                [1.0, 0.0, 0.0],
            ],
            cell=[20.0, 20.0, 20.0],
            pbc=True,
        )

        with self.assertWarnsRegex(
            RuntimeWarning,
            r'Multiple O atoms.*Selecting nearest O atom index 0',
        ):
            O_indices, shifts = find_bonded_oxygens(atoms, [2])

        self.assertEqual(O_indices, [0])
        np.testing.assert_allclose(shifts, [[0.0, 0.0, 0.0]])

    def test_write_OH_analysis(self):
        with open('fixtures/in.yaml') as f:
            param = yaml.safe_load(f)

        with tempfile.TemporaryDirectory() as tmp_dir:
            analysis_file = os.path.join(tmp_dir, 'OH_analysis.csv')
            in_file = os.path.join(tmp_dir, 'in.yaml')
            param['OH_analysis_out'] = analysis_file
            with open(in_file, 'w') as f:
                yaml.safe_dump(param, f)

            get_OH_dipoles(in_file)

            analysis = pd.read_csv(analysis_file)

        self.assertEqual(
            analysis.columns.tolist(),
            [
                'H_index', 'OH_type', 'OH_bond_length',
                'Hbond_lengths', 'Hbond_O_indices',
                'Wrapped_x', 'Wrapped_y', 'Wrapped_z',
            ],
        )
        self.assertEqual(analysis['H_index'].tolist(), [1, 3])
        self.assertEqual(analysis['OH_type'].tolist(), ['other', 'other'])
        np.testing.assert_allclose(analysis['OH_bond_length'], [1.0, 1.0])
        self.assertEqual(analysis['Hbond_lengths'].tolist(), ['[]', '[]'])
        self.assertEqual(analysis['Hbond_O_indices'].tolist(), ['[]', '[]'])
        np.testing.assert_allclose(
            analysis.loc[0, ['Wrapped_x', 'Wrapped_y', 'Wrapped_z']]
            .to_numpy(dtype=float),
            [2.0, 1.0, 1.0],
        )

    def test_OH_analysis_against_hand_fixture(self):
        fixture_dir = 'fixtures_b_OH_analysis'
        with open(os.path.join(fixture_dir, 'in.yaml')) as f:
            param = yaml.safe_load(f)

        with tempfile.TemporaryDirectory() as tmp_dir:
            analysis_file = os.path.join(tmp_dir, 'OH_analysis.csv')
            in_file = os.path.join(tmp_dir, 'in.yaml')
            param['OH_analysis_out'] = analysis_file
            with open(in_file, 'w') as f:
                yaml.safe_dump(param, f)

            get_OH_dipoles(in_file)
            actual = self.load_OH_analysis(analysis_file)

        expected = self.load_OH_analysis(
            os.path.join(fixture_dir, 'OH_analysis_hand.csv')
        )
        self.assertEqual(
            [
                (
                    row['H_index'],
                    row['OH_type'],
                    row['Hbond_O_indices'],
                )
                for row in actual
            ],
            [
                (
                    row['H_index'],
                    row['OH_type'],
                    row['Hbond_O_indices'],
                )
                for row in expected
            ],
        )
        np.testing.assert_allclose(
            [row['OH_bond_length'] for row in actual],
            [row['OH_bond_length'] for row in expected],
            atol=1e-5,
            rtol=0,
        )
        for actual_row, expected_row in zip(actual, expected):
            np.testing.assert_allclose(
                actual_row['Hbond_lengths'],
                expected_row['Hbond_lengths'],
                atol=1e-5,
                rtol=0,
            )
        np.testing.assert_allclose(
            [
                [row['Wrapped_x'], row['Wrapped_y'], row['Wrapped_z']]
                for row in actual
            ],
            [
                [row['Wrapped_x'], row['Wrapped_y'], row['Wrapped_z']]
                for row in expected
            ],
            atol=1e-5,
            rtol=0,
        )


if __name__ == '__main__':
    unittest.main()
