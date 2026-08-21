"""For partial IR analysis to be used later, make the trajectories of OH
dipole moments of selected Hs.

Input
    in.yaml, e.g.

        OH_dipole_H_indices : [700, 682]  # starts from 0
        # 700 : Si-OH, 682: Al-OH

        OH_dipole_H_labels : ['Si-OH', 'Al-OH']

        or

        OH_dipole_H_all : True  # for all H's

        dump_unfolded : '../a.traj/config.dump'
        atom_symbols : 'Si O H Al'

        charge : 'nominal'  # or 'born_isotropic' or 'born_full'
        born_file : 'xxx/sample_0_BORN'  # needed for 'born_isotropic' or 'born_full'
        born_poscar : 'xxx/POSCAR'       # needed for 'born_isotropic' or 'born_full'

        OH_dipoles_out : 'OH_dipoles.dat'   # optional, default 'OH_dipoles.dat'
        OH_analysis_out : 'OH_analysis.csv' # optional, default 'OH_analysis.csv'
        BORN_isotropic_out : 'BORN_isotropic.txt' # optional; null disables it

Output
    OH_dipoles.dat
        # step   (OH dipole components, |e| Ang)  H_700_x H_700_y H_700_z H_682_x H_682_y H_682_z  
        0    ...
        1    ...

    Each triplet is the OH bond dipole, which is the OH vector (pointing from
    the bonded O atom to the H atom) for the corresponding H index in
    OH_dipole_H_indices, multiplied by the charge to obtain the dipole.
    So the formula is H_charge * r(H-O), where H_charge is
    +1 for 'charge': 'nominal', and 1/3*(Z_11+Z_22+Z_33), the isotropic part
    for 'born_isotropic', and the 3x3 tensor Z_ij for 'born_full', where Z_ij is the
    Born effective charge tensor read by get_dipole_born.py.

    BORN_isotropic.txt: written for charge = 'born_isotropic'
        H_index  Born isotropic charge
        700      1.002

    OH_analysis.csv (can be read with pandas)
        H_index,OH_type,OH_bond_length,Hbond_lengths,Hbond_O_indices,Wrapped_x,Wrapped_y,Wrapped_z
    700,Si-OH,0.978,"[1.82]","[431]",12.34,15.67,8.90
    682,Al-OH,0.965,"[]","[]",11.20,14.50,7.80
    731,Si-OH,0.972,"[1.79,2.03]","[415,902]",13.10,16.20,9.40

    OH_type is one of ['Si-OH', 'Al-OH', 'Al..OH2', 'Al-(OH)-Si',
    'Al-(OH)-Al', 'Si-(OH)-Si', 'H2O', 'other'].
    The bond length is from the 0th configuration, which is supposed to be
    relaxed.
    For 'Al..OH2' or 'H2O', since the O is shared by two OH bonds, there are 
    stretching and antistretching modes. The H_index is the smaller numerical H atom index. The printed dipole moment
    is the sum, H1_charge * r(H1-O) + H2_charge * r(H2-O), where r(H1-O) is the
    vector from O to H1. This preserves the cross-correlation between the two bonds, allowing
    stretching and antistretching modes to have different IR intensities.

    The hydrogen_bond_length is the distance from H to a nearby O that is
    not directly bonded, for which the donor-O..acceptor-O distance is less
    than 3.5 Ang and the O-H..O angle is larger than 140 degrees. Otherwise
    mark the distance and O index as "[]". If many are found, write
    "[1.79,2.03]" and all acceptor O indices.
    For more info, please see the Supporting Info.
     H. Windeck, F. Berger, J. Sauer, Angew. Chem. Int. Ed. 2023, 62, e202303204. Spectroscopic Signatures of Internal Hydrogen Bonds of Brønsted-Acid Sites in the Zeolite H-MFI, https://doi.org/10.1002/anie.202303204
    The 3.5 Ang cutoff and 140 degrees is from J Phys Chem C Nanomater Interfaces. 2012 Nov 6;116(50):26247–26261. doi: 10.1021/jp302428b



Procedure
    It's similar to get_dipole_born.py, but instead of the total dipole moment, we
    track individual O-H bond vectors.

    1. Open dump_unfolded file (unwrapped/unfolded LAMMPS-style trajectory)
       via thermo_SiO2.io.read_sil.
    2. Select H atoms from OH_dipole_H_indices, or select every H atom in
       the first frame when OH_dipole_H_all is True. For each selected H,
       find the bonded O atom: the nearest O atom under the minimum image
       convention (MIC), with bond length < 1.5 Ang (typical O-H bond
       length is ~1.0 Ang). Print the found O index and the bond length.
       If no O atom is found within 1.5 Ang, raise an error and stop. If
       multiple O atoms are found, warn and select the nearest one.
    3. Since the dump file coordinates are unfolded (no periodic wrapping
       between frames) and computing the MIC every frame would be
       expensive, we instead compute, once, a constant Cartesian shift
       vector per O-H pair that maps the raw (unfolded) O position onto
       the same periodic image as the H atom, using the first frame's MIC
       vector as reference. This shift is then reused for every frame.
    4. Apply the requested charge to each O-H vector. For H2O and
       Al..OH2, automatically include both bonded H atoms and sum their
       dipoles before writing one triplet under the smaller H index.

    For OH_analysis.csv, OH_type is determined from the first-frame bond
    graph. H2O requires an O bonded only to two H atoms, with each H bonded
    only to that O. Other types use the Si and Al neighbors of the bonded
    O atom. Bonds are found with ASE's neighbor_list() using covalent
    radii multiplied by 1.5 (natural_cutoffs(mult=1.5)).
"""
import warnings

import numpy as np
import pandas as pd
import yaml
from ase.geometry import find_mic
from ase.io import read
from ase.neighborlist import natural_cutoffs, neighbor_list
from thermo_SiO2.IR_spectrum.dipole.get_dipole_born import read_born_charges
from thermo_SiO2.io import read_sil


OH_TYPES_BY_CATION_COUNTS = {
    (1, 0): 'Si-OH',
    (0, 1): 'Al-OH',
    (1, 1): 'Al-(OH)-Si',
    (0, 2): 'Al-(OH)-Al',
    (2, 0): 'Si-(OH)-Si',
}

BOND_CUTOFF_MULTIPLIER = 1.5
HYDROGEN_BOND_MAX_OO_DISTANCE = 3.5
HYDROGEN_BOND_MIN_ANGLE_DEG = 140.0


def find_bonded_oxygens(cfg0, H_indices, bond_cutoff=1.5):
    """For each H index, find the nearest O atom (via MIC) within
    bond_cutoff Ang. Returns:
        O_indices   : list of matched O atom index for each H
        shifts      : (n_pairs, 3) array; Cartesian shift to subtract from
                      the raw O position each frame so that
                      O_pos - shift is in the same periodic image as H_pos
    """
    symbols = np.array(cfg0.get_chemical_symbols())
    O_indices_all = np.where(symbols == 'O')[0]
    if O_indices_all.size == 0:
        raise ValueError('No O atoms found in the structure.')

    positions0 = cfg0.get_positions()
    cell = cfg0.cell
    pbc = cfg0.pbc

    O_indices = []
    shifts = []
    for h_idx in H_indices:
        if symbols[h_idx] != 'H':
            raise ValueError(
                f'Atom index {h_idx} in OH_dipole_H_indices is not an H '
                f'atom (found {symbols[h_idx]}).'
            )
        H_pos = positions0[h_idx]

        # vectors from H to every O atom, then take the minimum image
        D = positions0[O_indices_all] - H_pos
        mic_vectors, mic_dist = find_mic(D, cell=cell, pbc=pbc)

        matches = np.where(mic_dist < bond_cutoff)[0]
        if matches.size == 0:
            raise ValueError(
                f'No O atom found within {bond_cutoff} Ang of H atom '
                f'index {h_idx}.'
            )
        elif matches.size > 1:
            matched_O = O_indices_all[matches]
            matched_dist = mic_dist[matches]
            nearest = int(np.argmin(matched_dist))
            match = matches[nearest]
            warnings.warn(
                f'Multiple O atoms found within {bond_cutoff} Ang of H '
                f'atom index {h_idx}: O indices {matched_O.tolist()}, '
                f'distances {matched_dist.tolist()}. Selecting nearest '
                f'O atom index {int(matched_O[nearest])} with distance '
                f'{matched_dist[nearest]}.',
                RuntimeWarning,
                stacklevel=2,
            )
        else:
            match = matches[0]
        o_idx = int(O_indices_all[match])
        bond_len = mic_dist[match]
        mic_vec = mic_vectors[match]  # MIC vector from H to O (O - H)

        print(f'H index {h_idx}: bonded O index {o_idx}, '
              f'bond length {bond_len:.4f} Ang')

        # raw (unfolded) O - H difference at frame 0
        raw_diff = positions0[o_idx] - H_pos
        # constant lattice shift so that (O_raw - shift) - H = mic_vec
        shift = raw_diff - mic_vec

        O_indices.append(o_idx)
        shifts.append(shift)

    return O_indices, np.array(shifts)


def write_OH_dipoles(out_file, oh_dipoles, H_indices):
    """Write dipoles with shape (n_frames, n_groups, 3)."""
    header_cols = ' '.join(
        f'H_{h}_x H_{h}_y H_{h}_z' for h in H_indices
    )
    with open(out_file, 'w') as f:
        f.write(f'# step {header_cols}  (OH dipole components, |e| Ang)\n')
        for step, frame in enumerate(oh_dipoles):
            vals = ' '.join('{:.5f} {:.5f} {:.5f}'.format(*vec) for vec in frame)
            f.write(f'{step} {vals}\n')


def write_BORN_isotropic(out_file, H_indices, charge_tensors):
    """Write the isotropic Born charge used for each contributing H."""
    with open(out_file, 'w') as f:
        f.write('# H_index Born isotropic charge\n')
        for h_idx in H_indices:
            f.write(f'{h_idx} {charge_tensors[h_idx, 0, 0]:.5f}\n')


def classify_OH_types(cfg0, O_indices):
    """Classify OH groups from each O atom's first-frame Si/Al neighbors."""
    symbols = np.array(cfg0.get_chemical_symbols())
    neighbor_i, neighbor_j = neighbor_list(
        'ij',
        cfg0,
        cutoff=natural_cutoffs(cfg0, mult=BOND_CUTOFF_MULTIPLIER),
    )

    neighbors_by_atom = [set() for _ in range(len(cfg0))]
    for i, j in zip(neighbor_i, neighbor_j):
        neighbors_by_atom[int(i)].add(int(j))

    OH_types = []
    for o_idx in O_indices:
        neighbors = neighbors_by_atom[o_idx]
        neighbor_symbols = [symbols[j] for j in neighbors]
        hydrogen_neighbors = [
            j for j in neighbors if symbols[j] == 'H'
        ]
        cation_symbols = [
            symbol for symbol in neighbor_symbols
            if symbol in ('Si', 'Al')
        ]
        cation_counts = (
            cation_symbols.count('Si'),
            cation_symbols.count('Al'),
        )
        is_water = (
            len(neighbors) == 2
            and len(hydrogen_neighbors) == 2
            and all(
                neighbors_by_atom[h_idx] == {o_idx}
                for h_idx in hydrogen_neighbors
            )
        )
        is_Al_OH2 = (
            len(neighbors) == 3
            and len(hydrogen_neighbors) == 2
            and cation_counts == (0, 1)
        )
        if is_Al_OH2:
            OH_types.append('Al..OH2')
        elif is_water:
            OH_types.append('H2O')
        else:
            OH_types.append(
                OH_TYPES_BY_CATION_COUNTS.get(cation_counts, 'other')
            )

    return OH_types


def get_neighbors_by_atom(cfg0):
    """Return the first-frame covalent-neighbor graph."""
    neighbor_i, neighbor_j = neighbor_list(
        'ij',
        cfg0,
        cutoff=natural_cutoffs(cfg0, mult=BOND_CUTOFF_MULTIPLIER),
    )
    neighbors_by_atom = [set() for _ in range(len(cfg0))]
    for i, j in zip(neighbor_i, neighbor_j):
        neighbors_by_atom[int(i)].add(int(j))
    return neighbors_by_atom


def build_dipole_groups(cfg0, selected_H_indices):
    """Build output groups, combining H2O and Al..OH2 by shared O.

    Selecting either H of a two-H group selects the whole group. The smaller
    numerical H index is its stable output identifier.
    """
    selected_O_indices, selected_shifts = find_bonded_oxygens(
        cfg0, selected_H_indices
    )
    pair_by_H = {
        h_idx: (o_idx, shift)
        for h_idx, o_idx, shift in zip(
            selected_H_indices, selected_O_indices, selected_shifts
        )
    }
    neighbors_by_atom = get_neighbors_by_atom(cfg0)
    symbols = np.array(cfg0.get_chemical_symbols())

    groups = []
    seen_group_keys = set()
    for selected_h_idx, o_idx in zip(selected_H_indices, selected_O_indices):
        OH_type = classify_OH_types(cfg0, [o_idx])[0]
        if OH_type in ('H2O', 'Al..OH2'):
            H_indices = sorted(
                idx for idx in neighbors_by_atom[o_idx]
                if symbols[idx] == 'H'
            )
            if len(H_indices) != 2:
                raise ValueError(
                    f'{OH_type} at O index {o_idx} must have exactly two '
                    f'H neighbors; found {H_indices}.'
                )
            group_key = ('O', o_idx)
        else:
            H_indices = [selected_h_idx]
            group_key = ('H', selected_h_idx)

        if group_key in seen_group_keys:
            continue
        seen_group_keys.add(group_key)

        missing_H_indices = [
            h_idx for h_idx in H_indices if h_idx not in pair_by_H
        ]
        if missing_H_indices:
            added_O_indices, added_shifts = find_bonded_oxygens(
                cfg0, missing_H_indices
            )
            for h_idx, added_o_idx, shift in zip(
                missing_H_indices, added_O_indices, added_shifts
            ):
                if added_o_idx != o_idx:
                    raise ValueError(
                        f'H index {h_idx} was classified in the group at O '
                        f'index {o_idx}, but its nearest O is {added_o_idx}.'
                    )
                pair_by_H[h_idx] = (added_o_idx, shift)

        groups.append({
            'H_index': min(H_indices),
            'H_indices': H_indices,
            'O_index': o_idx,
            'OH_type': OH_type,
        })

    return groups, pair_by_H


def get_charge_tensors(param, cfg0):
    """Return one 3x3 charge tensor per atom for the requested model."""
    if 'charge' not in param:
        raise ValueError(
            "Missing required 'charge' setting. Set it to 'nominal', "
            "'born_isotropic', or 'born_full'."
        )
    charge = param['charge']
    valid_charges = ('nominal', 'born_isotropic', 'born_full')
    if charge not in valid_charges:
        raise ValueError(
            f'Invalid charge {charge!r}; expected one of {valid_charges}.'
        )
    if charge == 'nominal':
        return np.repeat(np.eye(3)[None, :, :], len(cfg0), axis=0)

    missing = [key for key in ('born_file', 'born_poscar') if not param.get(key)]
    if missing:
        raise ValueError(
            f"charge: {charge!r} requires {', '.join(missing)}."
        )

    cfg_born = read(param['born_poscar'])
    if len(cfg_born) != len(cfg0):
        raise ValueError(
            f'Number of atoms in born_poscar ({len(cfg_born)}) does not '
            f'match the trajectory ({len(cfg0)}).'
        )
    if not (cfg0.symbols == cfg_born.symbols).all():
        raise ValueError(
            'The atomic symbols/order of born_poscar and the trajectory differ.'
        )
    if not np.allclose(cfg0.cell.array, cfg_born.cell.array):
        raise ValueError('The cells of born_poscar and the trajectory differ.')

    dr = cfg0.positions - cfg_born.positions
    _, dr_dist = find_mic(dr, cell=cfg0.cell, pbc=cfg0.pbc)
    dr_max = float(np.max(dr_dist))
    if dr_max >= 5.0:
        raise ValueError(
            'The max displacement between born_poscar and the trajectory '
            f'is very large ({dr_max:.2f} Ang). Check structure identity '
            'and atom ordering.'
        )
    if dr_max > 3.0:
        warnings.warn(
            'The max displacement between born_poscar and the trajectory '
            f'is large ({dr_max:.2f} Ang). Check structure identity and '
            'atom ordering.',
            RuntimeWarning,
            stacklevel=2,
        )

    Z = read_born_charges(param['born_file'])
    if Z.shape != (len(cfg0), 3, 3):
        raise ValueError(
            f'Born tensors have shape {Z.shape}; expected '
            f'({len(cfg0)}, 3, 3).'
        )
    if charge == 'born_isotropic':
        isotropic = np.trace(Z, axis1=1, axis2=2) / 3.0
        Z = isotropic[:, None, None] * np.eye(3)[None, :, :]
    return Z


def get_OH_bond_lengths(cfg0, H_indices, O_indices, shifts):
    """Return first-frame MIC O-H bond lengths in Angstrom."""
    positions = cfg0.get_positions()
    bond_vectors = [
        positions[h_idx] - (positions[o_idx] - shift)
        for h_idx, o_idx, shift in zip(H_indices, O_indices, shifts)
    ]
    return np.linalg.norm(bond_vectors, axis=1)


def find_hydrogen_bonds(
    cfg0,
    H_indices,
    donor_O_indices,
    max_OO_distance=HYDROGEN_BOND_MAX_OO_DISTANCE,
    min_angle_deg=HYDROGEN_BOND_MIN_ANGLE_DEG,
):
    """Return qualifying (acceptor O index, H..O distance) pairs per H."""
    symbols = np.array(cfg0.get_chemical_symbols())
    all_O_indices = np.flatnonzero(symbols == 'O')
    positions = cfg0.get_positions()

    hydrogen_bonds = []
    for h_idx, donor_O_idx in zip(H_indices, donor_O_indices):
        acceptor_O_indices = all_O_indices[all_O_indices != donor_O_idx]
        donor_to_H, _ = find_mic(
            positions[h_idx] - positions[donor_O_idx],
            cell=cfg0.cell,
            pbc=cfg0.pbc,
        )
        donor_to_acceptors, OO_distances = find_mic(
            positions[acceptor_O_indices] - positions[donor_O_idx],
            cell=cfg0.cell,
            pbc=cfg0.pbc,
        )

        H_to_donor = -donor_to_H
        H_to_acceptors = donor_to_acceptors - donor_to_H
        H_to_acceptor_distances = np.linalg.norm(H_to_acceptors, axis=1)
        cos_angles = np.einsum(
            'ij,j->i', H_to_acceptors, H_to_donor
        ) / (H_to_acceptor_distances * np.linalg.norm(H_to_donor))
        angles = np.degrees(np.arccos(np.clip(cos_angles, -1.0, 1.0)))

        matches = [
            (int(acceptor_O_idx), float(HO_distance))
            for acceptor_O_idx, HO_distance, OO_distance, angle in zip(
                acceptor_O_indices,
                H_to_acceptor_distances,
                OO_distances,
                angles,
            )
            if OO_distance < max_OO_distance and angle > min_angle_deg
        ]
        matches.sort(key=lambda match: (match[1], match[0]))
        hydrogen_bonds.append(matches)

    return hydrogen_bonds


def write_OH_analysis(
    out_file,
    H_indices,
    OH_types,
    bond_lengths,
    hydrogen_bonds,
    wrapped_H_positions,
):
    """Write first-frame OH and hydrogen-bond analysis as CSV."""
    rows = []
    for h_idx, OH_type, bond_length, H_bonds, H_position in zip(
        H_indices,
        OH_types,
        bond_lengths,
        hydrogen_bonds,
        wrapped_H_positions,
    ):
        x, y, z = H_position
        rows.append({
            'H_index': h_idx,
            'OH_type': OH_type,
            'OH_bond_length': bond_length,
            'Hbond_lengths': [
                round(distance, 5) for _, distance in H_bonds
            ],
            'Hbond_O_indices': [O_idx for O_idx, _ in H_bonds],
            'Wrapped_x': x,
            'Wrapped_y': y,
            'Wrapped_z': z,
        })

    pd.DataFrame(rows).to_csv(out_file, index=False, float_format='%.5f')


def select_H_indices(param, cfg0):
    """Select explicit H indices or all H indices from the first frame."""
    has_explicit_indices = 'OH_dipole_H_indices' in param
    select_all = param.get('OH_dipole_H_all', False)

    if has_explicit_indices and select_all:
        raise ValueError(
            'Set either OH_dipole_H_indices or OH_dipole_H_all: True, '
            'not both.'
        )
    if select_all:
        symbols = np.array(cfg0.get_chemical_symbols())
        return np.flatnonzero(symbols == 'H').tolist()
    if has_explicit_indices:
        return param['OH_dipole_H_indices']
    raise ValueError(
        'Set OH_dipole_H_indices or OH_dipole_H_all: True in the input.'
    )


def compute_frame_dipoles(positions, groups, pair_by_H, charge_tensors):
    """Compute all selected OH/group dipoles for one trajectory frame."""
    frame_dipoles = []
    for group in groups:
        dipole = np.zeros(3)
        for h_idx in group['H_indices']:
            o_idx, shift = pair_by_H[h_idx]
            oh_vector = positions[h_idx] - (positions[o_idx] - shift)
            dipole += charge_tensors[h_idx] @ oh_vector
        frame_dipoles.append(dipole)
    return frame_dipoles


def get_OH_dipoles(in_file='in.yaml'):
    with open(in_file, 'r') as stream:
        try:
            param = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
            raise

    if 'charge' not in param:
        raise ValueError(
            "Missing required 'charge' setting. Set it to 'nominal', "
            "'born_isotropic', or 'born_full'."
        )

    cfgs = read_sil(
        param['dump_unfolded'], atom_symbols=param['atom_symbols']
    )
    print('cfgs were read.')

    OH_dipoles_out = param.get('OH_dipoles_out', 'OH_dipoles.dat')
    OH_analysis_out = param.get('OH_analysis_out', 'OH_analysis.csv')
    BORN_isotropic_out = param.get(
        'BORN_isotropic_out', 'BORN_isotropic.txt'
    )

    cfg0 = cfgs[0]
    selected_H_indices = select_H_indices(param, cfg0)
    groups, pair_by_H = build_dipole_groups(cfg0, selected_H_indices)
    charge_tensors = get_charge_tensors(param, cfg0)
    output_H_indices = [group['H_index'] for group in groups]

    if (
        param['charge'] == 'born_isotropic'
        and BORN_isotropic_out is not None
    ):
        contributing_H_indices = [
            h_idx
            for group in groups
            for h_idx in group['H_indices']
        ]
        write_BORN_isotropic(
            BORN_isotropic_out,
            contributing_H_indices,
            charge_tensors,
        )

    if OH_analysis_out is not None:
        representative_O_indices = [group['O_index'] for group in groups]
        representative_shifts = [
            pair_by_H[group['H_index']][1] for group in groups
        ]
        bond_lengths = get_OH_bond_lengths(
            cfg0,
            output_H_indices,
            representative_O_indices,
            representative_shifts,
        )
        hydrogen_bonds = find_hydrogen_bonds(
            cfg0, output_H_indices, representative_O_indices
        )
        write_OH_analysis(
            OH_analysis_out,
            output_H_indices,
            [group['OH_type'] for group in groups],
            bond_lengths,
            hydrogen_bonds,
            cfg0.get_positions(wrap=True)[output_H_indices],
        )

    print(f'total {len(cfgs)} cfgs')
    oh_dipoles = []
    for i, atoms in enumerate(cfgs):
        if i % 1000 == 0:
            print(f'{i} / {len(cfgs)} cfgs')
        positions = atoms.get_positions()

        frame_dipoles = compute_frame_dipoles(
            positions, groups, pair_by_H, charge_tensors
        )
        oh_dipoles.append(frame_dipoles)

    oh_dipoles = np.asarray(oh_dipoles)

    if OH_dipoles_out is not None:
        write_OH_dipoles(OH_dipoles_out, oh_dipoles, output_H_indices)

    return oh_dipoles


if __name__ == '__main__':
    get_OH_dipoles()
