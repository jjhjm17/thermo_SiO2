"""For partial IR analysis to be used later, make the trajectories of OH
vectors of selected Hs.

Input
    in.yaml, e.g.

        OH_vector_H_indices : [700, 682]  # starts from 0
        # 700 : Si-OH, 682: Al-OH
        OH_vector_H_labels : ['Si-OH', 'Al-OH']

        or

        OH_vector_H_all : True  # for all H's

        dump_unfolded : '../a.traj/config.dump'
        atom_symbols : 'Si O H Al'
        OH_vectors_out : 'OH_vectors.dat'   # optional, default 'OH_vectors.dat'
        OH_analysis_out : 'OH_analysis.dat' # optional, default 'OH_analysis.dat'

Output
    OH_vectors.dat
        # step  H_700_x H_700_y H_700_z   H_682_x H_682_y H_682_z
        0    ...
        1    ...

    Each triplet is the OH bond vector (pointing from the bonded O atom to
    the H atom) for the corresponding H index in OH_vector_H_indices.

    OH_analysis.dat
       # H_index OH_type OH_bond_length (Ang) Wrapped_x Wrapped_y Wrapped_z (Ang)
       700 Si-OH 0.9xx ... ... ...
       682 Al-OH 0.9xx ... ... ...

    OH_type is one of ['Si-OH', 'Al-OH', 'Al..OH2', 'Al-(OH)-Si',
    'Al-(OH)-Al', 'Si-(OH)-Si', 'H2O', 'other'].
    The bond length is from the 0th configuration, which is supposed to be
    relaxed.

Procedure
    It's similar to get_dipole_born.py, but instead of Born charges, we
    track individual O-H bond vectors.

    1. Open dump_unfolded file (unwrapped/unfolded LAMMPS-style trajectory)
       via thermo_SiO2.io.read_sil.
    2. Select H atoms from OH_vector_H_indices, or select every H atom in
       the first frame when OH_vector_H_all is True. For each selected H,
       find the bonded O atom: the nearest O atom under the minimum image
       convention (MIC), with bond length < 1.5 Ang (typical O-H bond
       length is ~1.0 Ang). Print the found O index and the bond length.
       If zero or more than one O atom is found within 1.5 Ang, raise an
       error and stop.
    3. Since the dump file coordinates are unfolded (no periodic wrapping
       between frames) and computing the MIC every frame would be
       expensive, we instead compute, once, a constant Cartesian shift
       vector per O-H pair that maps the raw (unfolded) O position onto
       the same periodic image as the H atom, using the first frame's MIC
       vector as reference. This shift is then reused for every frame.
    4. Loop over all structures in the dump file. For each H atom in
       OH_vector_H_indices, apply the pair's stored shift to the O
       position ("wrap" it), then compute OH vector = H_pos - O_wrapped,
       and write out its components.

    For OH_analysis.dat, OH_type is determined from the first-frame bond
    graph. H2O requires an O bonded only to two H atoms, with each H bonded
    only to that O. Other types use the Si and Al neighbors of the bonded
    O atom. Bonds are found with ASE's neighbor_list() using covalent
    radii multiplied by 1.5 (natural_cutoffs(mult=1.5)).
"""
import numpy as np
import yaml
from ase.geometry import find_mic
from ase.neighborlist import natural_cutoffs, neighbor_list
from thermo_SiO2.io import read_sil


OH_TYPES_BY_CATION_COUNTS = {
    (1, 0): 'Si-OH',
    (0, 1): 'Al-OH',
    (1, 1): 'Al-(OH)-Si',
    (0, 2): 'Al-(OH)-Al',
    (2, 0): 'Si-(OH)-Si',
}

BOND_CUTOFF_MULTIPLIER = 1.5


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
                f'Atom index {h_idx} in OH_vector_H_indices is not an H '
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
            raise ValueError(
                f'Multiple O atoms found within {bond_cutoff} Ang of H '
                f'atom index {h_idx}: O indices {matched_O.tolist()}, '
                f'distances {matched_dist.tolist()}.'
            )

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


def write_OH_vectors(out_file, oh_vectors, H_indices):
    """oh_vectors: shape (n_frames, n_pairs, 3)"""
    header_cols = ' '.join(
        f'H_{h}_x H_{h}_y H_{h}_z' for h in H_indices
    )
    with open(out_file, 'w') as f:
        f.write(f'# step {header_cols}  (OH vector components, Ang)\n')
        for step, frame in enumerate(oh_vectors):
            vals = ' '.join('{:.5f} {:.5f} {:.5f}'.format(*vec) for vec in frame)
            f.write(f'{step} {vals}\n')


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


def get_OH_bond_lengths(cfg0, H_indices, O_indices, shifts):
    """Return first-frame MIC O-H bond lengths in Angstrom."""
    positions = cfg0.get_positions()
    bond_vectors = [
        positions[h_idx] - (positions[o_idx] - shift)
        for h_idx, o_idx, shift in zip(H_indices, O_indices, shifts)
    ]
    return np.linalg.norm(bond_vectors, axis=1)


def write_OH_analysis(
    out_file, H_indices, OH_types, bond_lengths, wrapped_H_positions
):
    """Write first-frame OH topology, lengths, and wrapped H positions."""
    with open(out_file, 'w') as f:
        f.write(
            '# H_index OH_type OH_bond_length (Ang) '
            'Wrapped_x Wrapped_y Wrapped_z (Ang)\n'
        )
        for h_idx, OH_type, bond_length, H_position in zip(
            H_indices, OH_types, bond_lengths, wrapped_H_positions
        ):
            x, y, z = H_position
            f.write(
                f'{h_idx} {OH_type} {bond_length:.5f} '
                f'{x:.5f} {y:.5f} {z:.5f}\n'
            )


def select_H_indices(param, cfg0):
    """Select explicit H indices or all H indices from the first frame."""
    has_explicit_indices = 'OH_vector_H_indices' in param
    select_all = param.get('OH_vector_H_all', False)

    if has_explicit_indices and select_all:
        raise ValueError(
            'Set either OH_vector_H_indices or OH_vector_H_all: True, '
            'not both.'
        )
    if select_all:
        symbols = np.array(cfg0.get_chemical_symbols())
        return np.flatnonzero(symbols == 'H').tolist()
    if has_explicit_indices:
        return param['OH_vector_H_indices']
    raise ValueError(
        'Set OH_vector_H_indices or OH_vector_H_all: True in the input.'
    )


def get_OH_vectors(in_file='in.yaml'):
    with open(in_file, 'r') as stream:
        try:
            param = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
            raise

    cfgs = read_sil(param['dump_unfolded'],
                     atom_symbols=param['atom_symbols'])
    print('cfgs were read.')

    OH_vectors_out = param.get('OH_vectors_out', 'OH_vectors.dat')
    OH_analysis_out = param.get('OH_analysis_out', 'OH_analysis.dat')

    cfg0 = cfgs[0]
    H_indices = select_H_indices(param, cfg0)
    O_indices, shifts = find_bonded_oxygens(cfg0, H_indices)

    if OH_analysis_out is not None:
        OH_types = classify_OH_types(cfg0, O_indices)
        bond_lengths = get_OH_bond_lengths(
            cfg0, H_indices, O_indices, shifts
        )
        write_OH_analysis(
            OH_analysis_out,
            H_indices,
            OH_types,
            bond_lengths,
            cfg0.get_positions(wrap=True)[H_indices],
        )

    print(f'total {len(cfgs)} cfgs')
    oh_vectors = []
    for i, atoms in enumerate(cfgs):
        if i % 1000 == 0:
            print(f'{i} / {len(cfgs)} cfgs')
        positions = atoms.get_positions()

        frame_vectors = []
        for h_idx, o_idx, shift in zip(H_indices, O_indices, shifts):
            H_pos = positions[h_idx]
            O_wrap = positions[o_idx] - shift
            oh_vec = H_pos - O_wrap
            frame_vectors.append(oh_vec)

        oh_vectors.append(frame_vectors)

    oh_vectors = np.array(oh_vectors)  # (n_frames, n_pairs, 3)

    if OH_vectors_out is not None:
        write_OH_vectors(OH_vectors_out, oh_vectors, H_indices)

    return oh_vectors


if __name__ == '__main__':
    get_OH_vectors()
