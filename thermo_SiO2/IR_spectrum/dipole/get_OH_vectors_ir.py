"""
This script obtains IR spectrum of each OH groups.
We use formal charge of H: +1 to get the dipole moment from the OH vector,
which depends on the difference of OH positions.
Input
    in.yaml:
        OH_vector_H_indices : [700, 682]  # starts from 0
        or
        OH_vector_H_all : True  # read all H indices from OH_vectors.dat
        OH_output_plot = "IR_spectrum_OH.pdf"
        OH_output_data = "IR_spectrum_OH.dat"
        OH_analysis_out = "OH_analysis.dat"
        OH_normalize : False  # optional, default False; True normalizes
            each partial spectrum to its own max = 1 (as ir.py does),
            which makes different O-H bonds' peak heights NOT comparable

        # group
        OH_output_group_plot = "IR_spectrum_OH_group.pdf"
        OH_output_group_data = "IR_spectrum_OH_group.dat"
        OH_analysis_group_data : 'OH_analysis_group.dat'

    OH_vectors.dat
Output
    OH_output_plot: a plot of all partial IR's of OH vectors,
        with legens indicating the index of H atoms,
        similar to output_plot of ir.py.
    OH_output_data: data file, similar to the output_data of ir.py,
        but for many intensities for each H in OH_vector_H_indices.

    group
    OH_output_group_plot: a plot of partial IRs grouped by OH types.
        For each group, the sum of partial IRs' are plotted and labeled.
    OH_output_group_data: data file with one summed intensity column per
        OH type.

Details 
    After running get_OH_vectors.py, we have OH_vectors.dat
        ) head OH_vectors.dat
        # step H_700_x H_700_y H_700_z H_682_x H_682_y H_682_z  (OH vector components, Ang)
        0 0.57550 0.73990 0.22780 -0.93727 0.17790 -0.08117
        1 0.57575 0.75050 0.23230 -0.93041 0.18010 -0.08376
        ...

Procedure
    We use the same dipole-autocorrelation / FFT method as ir.py
    (ir.compute_ir_from_dipole_components), applied to each O-H bond
    separately, instead of to the whole-cell dipole moment.

    1. Read in.yaml for OH_vector_H_indices, or read all H indices in
       column order from the OH_vectors.dat header when OH_vector_H_all
       is True. Also read OH_vectors_out (input file written by
       get_OH_vectors.py, default 'OH_vectors.dat'), dt_fs, and the output
       filenames OH_output_plot / OH_output_data.
    2. Load OH_vectors.dat. For each H index, take its O-H vector
       (x, y, z) as the "dipole" components, scaled by the formal charge
       of H (+1 |e|), and run it through
       ir.compute_ir_from_dipole_components to get a partial IR
       spectrum.
    3. Save all partial spectra (sharing one frequency axis) to
       OH_output_data, one intensity column per H index.
       By default the spectra are NOT individually normalized (unlike
       ir.py's own output), so that the relative peak heights of
       different O-H bonds' spectra can be directly compared; set
       OH_normalize: True in in.yaml to instead normalize each partial
       spectrum to its own max = 1, as ir.py does.
    4. Plot all partial spectra on one figure, one line per H index,
       with a legend giving the H atom index, and save to
       OH_output_plot.
    5. Read each H atom's OH type from OH_analysis_out and sum the partial
       spectra within each type. Save the grouped plot to
       OH_output_group_plot and the grouped columns to
       OH_output_group_data. Set either output to null to disable it.
"""
import re
import numpy as np
import matplotlib.pyplot as plt
import yaml
from thermo_SiO2.IR_spectrum.dipole.ir import compute_ir_from_dipole_components

# Formal charge assigned to H when converting an O-H bond vector into a
# "dipole moment" for the FFT/autocorrelation method (see module
# docstring). With OH_normalize: False (the default here), this constant
# scales the absolute intensity of every partial spectrum by the same
# factor, so it does not affect their *relative* heights; it matters
# only if you want the absolute intensity to correspond to a real
# dipole moment (e.g. to compare against ir.py's whole-cell spectrum,
# also with normalize=False).
H_FORMAL_CHARGE = 1.0


def read_OH_vector_H_indices(OH_vectors_file):
    """Read and validate H-index xyz triplets from OH_vectors.dat."""
    with open(OH_vectors_file) as f:
        header = f.readline().strip()

    columns = re.findall(r'\bH_(\d+)_([xyz])\b', header)
    if not columns or len(columns) % 3 != 0:
        raise ValueError(
            f'Could not read complete H-index xyz triplets from the '
            f'header of {OH_vectors_file!r}.'
        )

    H_indices = []
    for start in range(0, len(columns), 3):
        triplet = columns[start:start + 3]
        indices = [int(h_idx) for h_idx, _ in triplet]
        axes = [axis for _, axis in triplet]
        if len(set(indices)) != 1 or axes != ['x', 'y', 'z']:
            raise ValueError(
                f'Invalid OH-vector column triplet {triplet!r} in the '
                f'header of {OH_vectors_file!r}.'
            )
        H_indices.append(indices[0])

    if len(set(H_indices)) != len(H_indices):
        raise ValueError(
            f'Duplicate H indices in the header of {OH_vectors_file!r}: '
            f'{H_indices}.'
        )
    return H_indices


def select_H_indices(param, OH_vectors_file):
    """Select explicit H indices or all indices from the vector header."""
    has_explicit_indices = 'OH_vector_H_indices' in param
    select_all = param.get('OH_vector_H_all', False)

    if has_explicit_indices and select_all:
        raise ValueError(
            'Set either OH_vector_H_indices or OH_vector_H_all: True, '
            'not both.'
        )
    if select_all:
        return read_OH_vector_H_indices(OH_vectors_file)
    if has_explicit_indices:
        return param['OH_vector_H_indices']
    raise ValueError(
        'Set OH_vector_H_indices or OH_vector_H_all: True in the input.'
    )


def load_OH_vectors(OH_vectors_file, H_indices):
    """Load OH_vectors.dat (as written by get_OH_vectors.py) and return
    a dict mapping each H index to its (N_frames, 3) array of O-H bond
    vector components (Ang).
    """
    data = np.loadtxt(OH_vectors_file)
    n_pairs = (data.shape[1] - 1) // 3
    if n_pairs != len(H_indices):
        raise ValueError(
            f'OH_vectors_file {OH_vectors_file!r} has {n_pairs} O-H '
            f'pairs, but OH_vector_H_indices in in.yaml lists '
            f'{len(H_indices)} indices: {H_indices}.'
        )

    oh_vectors = {}
    for i, h_idx in enumerate(H_indices):
        cols = slice(1 + 3 * i, 1 + 3 * i + 3)
        oh_vectors[h_idx] = data[:, cols]

    return oh_vectors


def load_OH_types(OH_analysis_file, H_indices):
    """Load the OH type for every selected H from OH_analysis.dat."""
    OH_types_by_H = {}
    with open(OH_analysis_file) as f:
        for line_number, line in enumerate(f, start=1):
            if line.startswith('#') or not line.strip():
                continue
            columns = line.split()
            if len(columns) < 2:
                raise ValueError(
                    f'Invalid row {line_number} in '
                    f'{OH_analysis_file!r}: {line.rstrip()!r}.'
                )
            try:
                h_idx = int(columns[0])
            except ValueError as exc:
                raise ValueError(
                    f'Invalid H index on row {line_number} in '
                    f'{OH_analysis_file!r}: {columns[0]!r}.'
                ) from exc
            if h_idx in OH_types_by_H:
                raise ValueError(
                    f'Duplicate H index {h_idx} in '
                    f'{OH_analysis_file!r}.'
                )
            OH_types_by_H[h_idx] = columns[1]

    missing = [h_idx for h_idx in H_indices if h_idx not in OH_types_by_H]
    if missing:
        raise ValueError(
            f'OH analysis file {OH_analysis_file!r} has no OH type for '
            f'H indices {missing}.'
        )
    return {h_idx: OH_types_by_H[h_idx] for h_idx in H_indices}


def sum_ir_by_OH_type(ir_by_H, OH_types_by_H):
    """Sum partial IR spectra by OH type, preserving first-seen order."""
    ir_by_OH_type = {}
    for h_idx, ir in ir_by_H.items():
        OH_type = OH_types_by_H[h_idx]
        if OH_type not in ir_by_OH_type:
            ir_by_OH_type[OH_type] = np.zeros_like(ir)
        ir_by_OH_type[OH_type] += ir
    return ir_by_OH_type


def write_OH_ir_data(out_file, freq_cm, ir_by_H):
    """ir_by_H: dict H_index -> 1D intensity array (same length as
    freq_cm, and in the same order as the dict iterates)."""
    H_indices = list(ir_by_H.keys())
    header_cols = ' '.join(f'H_{h}' for h in H_indices)
    columns = np.column_stack([freq_cm] + [ir_by_H[h] for h in H_indices])
    np.savetxt(out_file, columns,
               header=f'Frequency(cm^-1) {header_cols}', fmt='%.6g')


def write_OH_group_ir_data(out_file, freq_cm, ir_by_OH_type):
    """Write one summed intensity column for each OH type."""
    OH_types = list(ir_by_OH_type)
    columns = np.column_stack(
        [freq_cm] + [ir_by_OH_type[OH_type] for OH_type in OH_types]
    )
    header_cols = ' '.join(OH_types)
    np.savetxt(
        out_file,
        columns,
        header=f'Frequency(cm^-1) {header_cols}',
        fmt='%.6g',
    )


def plot_OH_ir(out_file, freq_cm, ir_by_H, OH_vector_H_labels=False, freq_range=None, normalize=False):
    plt.figure(figsize=(7, 5))
    for label_idx, (h_idx, ir) in enumerate(ir_by_H.items()):
        if OH_vector_H_labels is not False:
            label = OH_vector_H_labels[label_idx]
        else:
            label=f'H {h_idx}'
        plt.plot(freq_cm, ir, linewidth=1.5, label=label)
    if freq_range is not None:
        plt.xlim(freq_range)
    if normalize:
        plt.ylim([0, 1.05])
    else:
        plt.ylim(bottom=0)  # let matplotlib pick the top, shared across curves
    plt.xlabel("Wavenumber (cm$^{-1}$)")
    plt.ylabel("Intensity (AU)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_file, dpi=300)


def plot_OH_group_ir(out_file, freq_cm, ir_by_OH_type, freq_range=None):
    """Plot the summed partial IR spectrum for each OH type."""
    plt.figure(figsize=(7, 5))
    for OH_type, ir in ir_by_OH_type.items():
        plt.plot(freq_cm, ir, linewidth=1.5, label=OH_type)
    if freq_range is not None:
        plt.xlim(freq_range)
    plt.ylim(bottom=0)
    plt.xlabel("Wavenumber (cm$^{-1}$)")
    plt.ylabel("Intensity (AU)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_file, dpi=300)


def get_OH_vectors_ir(in_file='in.yaml'):
    with open(in_file, 'r') as stream:
        try:
            param = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
            raise

    OH_vectors_file = param.get('OH_vectors_out', 'OH_vectors.dat')
    H_indices = select_H_indices(param, OH_vectors_file)
    OH_output_plot = param.get('OH_output_plot', 'IR_spectrum_OH.pdf')
    OH_output_group_plot = param.get(
        'OH_output_group_plot', 'IR_spectrum_OH_group.pdf'
    )
    OH_output_group_data = param.get(
        'OH_output_group_data', 'IR_spectrum_OH_group.dat'
    )
    OH_output_data = param.get('OH_output_data', 'IR_spectrum_OH.dat')
    OH_analysis_file = param.get('OH_analysis_out', 'OH_analysis.dat')
    dt_fs = param['dt_fs']

    use_window = param.get('use_window', True)
    use_gradient = param.get('use_gradient', True)  # d mu / d t
    cut_autocorr = param.get('cut_autocorr', 1000)  # fs  # see MACE4IR paper
    freq_range = param.get('freq_range', None)
    # False by default: keep each partial spectrum's raw (unscaled)
    # intensity, so the different O-H bonds' peak heights stay directly
    # comparable to each other. Set OH_normalize: True in in.yaml to
    # instead rescale each spectrum to its own max = 1 (as ir.py does).
    normalize = param.get('OH_normalize', False)
    OH_vector_H_labels = param.get('OH_vector_H_labels', False)

    oh_vectors = load_OH_vectors(OH_vectors_file, H_indices)

    freq_cm = None
    ir_by_H = {}
    for h_idx in H_indices:
        vec = oh_vectors[h_idx]  # (N_frames, 3)
        # [:-1]: keep the length even, matching ir.py's convention
        # (a slightly faster FFT); harmless if the length was already
        # even.
        vec = vec[:-1]
        mu_x = H_FORMAL_CHARGE * vec[:, 0]
        mu_y = H_FORMAL_CHARGE * vec[:, 1]
        mu_z = H_FORMAL_CHARGE * vec[:, 2]

        freq_cm_h, ir_h = compute_ir_from_dipole_components(
            [mu_x, mu_y, mu_z], dt_fs,
            use_window=use_window, use_gradient=use_gradient,
            cut_autocorr=cut_autocorr, normalize=normalize,
        )
        if freq_cm is None:
            freq_cm = freq_cm_h
        ir_by_H[h_idx] = ir_h

    write_OH_ir_data(OH_output_data, freq_cm, ir_by_H)
    plot_OH_ir(OH_output_plot, freq_cm, ir_by_H, OH_vector_H_labels,
               freq_range=freq_range,
               normalize=normalize)

    if OH_output_group_plot is not None or OH_output_group_data is not None:
        OH_types_by_H = load_OH_types(OH_analysis_file, H_indices)
        ir_by_OH_type = sum_ir_by_OH_type(ir_by_H, OH_types_by_H)
        if OH_output_group_data is not None:
            write_OH_group_ir_data(
                OH_output_group_data, freq_cm, ir_by_OH_type
            )
        if OH_output_group_plot is not None:
            plot_OH_group_ir(
                OH_output_group_plot,
                freq_cm,
                ir_by_OH_type,
                freq_range=freq_range,
            )

    return freq_cm, ir_by_H


if __name__ == '__main__':
    get_OH_vectors_ir()
