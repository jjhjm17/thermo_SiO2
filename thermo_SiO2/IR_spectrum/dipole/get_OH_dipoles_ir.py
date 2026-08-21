"""
This script obtains IR spectrum of each OH groups.

Input
    in.yaml:
        OH_dipole_H_indices : [700, 682]  # starts from 0
        or
        OH_dipole_H_all : True  # read all H indices from OH_dipoles.dat
        OH_output_plot : "IR_spectrum_OH.pdf"
        OH_output_data : "IR_spectrum_OH.dat"
        OH_analysis_out : "OH_analysis.csv"
        OH_normalize : False  # optional, default False; True normalizes
            each partial spectrum to its own max = 1 (as ir.py does),
            which makes different O-H bonds' peak heights NOT comparable

        # group
        OH_output_group_plot : "IR_spectrum_OH_group.pdf"
        OH_output_group_data : "IR_spectrum_OH_group.dat"
        OH_analysis_group_data : 'OH_analysis_group.dat'
        OH_freq_vs_bond_length_plot : 'OH_freq_vs_bond_length.pdf'
        OH_freq_data : 'OH_freq.dat'

    OH_dipoles.dat

Output
    OH_output_plot: a plot of all partial IR's of OH dipoles,
        with legens indicating the index of H atoms,
        similar to output_plot of ir.py.
    OH_output_data: data file, similar to the output_data of ir.py,
        but for many intensities for each H in OH_dipole_H_indices.

    group
    OH_output_group_plot: a plot of partial IRs grouped by OH types.
        For each group, the sum of partial IRs' are plotted and labeled.
    OH_output_group_data: data file with one summed intensity column per
        OH type.

    OH_freq_vs_bond_length_plot, OH_freq_data :
        We expect a linear relationship between OH frequency and OH bond length.

        The format:
            # H_index  OH_freq (cm-1)  OH_freq_stdev (cm-1)
            682        3700.xx          100.xx

        The OH stretching frequency is obtained as the average of frequencies over
        2000 cm-1, weighted by the IR intensity. The standard deviation is also calculated.

        For the plotting, first group all data by OH_type in OH_analysis_out file, and also
        group by 'free' or 'H-bonded'. Here, 'H-bonded' means there is at least one atom in
        Hbond_O_indices, and 'free' means none. Use filled symbols for 'H-bonded' and empty
        symbols for 'free' OH groups.

        Add a linear regression line and write the R2 value and the equation in the plot.
        The regression is done globally and not per OH type.

        See Fig. S5 Windeck et al, 2023
        https://doi.org/10.1002/anie.202303204
        Angewandte Chemie International Edition 62, Issue 25 e202303204

Details 
    After running get_OH_dipoles.py, we have OH_dipoles.dat
        ) head OH_dipoles.dat
        # step   (OH dipole components, |e| Ang)  H_700_x H_700_y H_700_z H_682_x H_682_y H_682_z  
        0 0.57550 0.73990 0.22780 -0.93727 0.17790 -0.08117
        1 0.57575 0.75050 0.23230 -0.93041 0.18010 -0.08376
        ...

Procedure
    We use the same dipole-autocorrelation / FFT method as ir.py
    (ir.compute_ir_from_dipole_components), applied to each stored OH or two-H molecular group, instead of to the whole-cell dipole moment.

    1. Read in.yaml for OH_dipole_H_indices, or read all H indices in
       column order from the OH_dipoles.dat header when OH_dipole_H_all
       is True. Also read OH_dipoles_out (input file written by
       get_OH_dipoles.py, default 'OH_dipoles.dat'), dt_fs, and the output
       filenames OH_output_plot / OH_output_data.
    2. Load OH_dipoles.dat. For each representative H index, take its stored OH or molecular dipole
       (x, y, z) as the "dipole" components. and run it through
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
import ast
import re
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import yaml
from thermo_SiO2.IR_spectrum.dipole.ir import compute_ir_from_dipole_components

OH_STRETCH_MIN_FREQ_CM = 2000.0


def read_OH_dipole_H_indices(OH_dipoles_file):
    """Read and validate H-index xyz triplets from OH_dipoles.dat."""
    with open(OH_dipoles_file) as f:
        header = f.readline().strip()

    columns = re.findall(r'\bH_(\d+)_([xyz])\b', header)
    if not columns or len(columns) % 3 != 0:
        raise ValueError(
            f'Could not read complete H-index xyz triplets from the '
            f'header of {OH_dipoles_file!r}.'
        )

    H_indices = []
    for start in range(0, len(columns), 3):
        triplet = columns[start:start + 3]
        indices = [int(h_idx) for h_idx, _ in triplet]
        axes = [axis for _, axis in triplet]
        if len(set(indices)) != 1 or axes != ['x', 'y', 'z']:
            raise ValueError(
                f'Invalid OH-dipole column triplet {triplet!r} in the '
                f'header of {OH_dipoles_file!r}.'
            )
        H_indices.append(indices[0])

    if len(set(H_indices)) != len(H_indices):
        raise ValueError(
            f'Duplicate H indices in the header of {OH_dipoles_file!r}: '
            f'{H_indices}.'
        )
    return H_indices


def select_H_indices(param, OH_dipoles_file):
    """Validate selection settings and return the serialized group indices.

    The file header is authoritative because H2O and Al..OH2 selections are
    deduplicated and serialized under the smaller of their two H indices.
    """
    has_explicit_indices = 'OH_dipole_H_indices' in param
    select_all = param.get('OH_dipole_H_all', False)

    if has_explicit_indices and select_all:
        raise ValueError(
            'Set either OH_dipole_H_indices or OH_dipole_H_all: True, '
            'not both.'
        )
    if not has_explicit_indices and not select_all:
        raise ValueError(
            'Set OH_dipole_H_indices or OH_dipole_H_all: True in the input.'
        )
    return read_OH_dipole_H_indices(OH_dipoles_file)


def load_OH_dipoles(OH_dipoles_file, H_indices):
    """Load OH_dipoles.dat (as written by get_OH_dipoles.py) and return
    a dict mapping each H index to its (N_frames, 3) array of OH or molecular-group
    dipole components (|e| Ang).
    """
    data = np.loadtxt(OH_dipoles_file)
    n_pairs = (data.shape[1] - 1) // 3
    if n_pairs != len(H_indices):
        raise ValueError(
            f'OH_dipoles_file {OH_dipoles_file!r} has {n_pairs} O-H '
            f'pairs, but OH_dipole_H_indices in in.yaml lists '
            f'{len(H_indices)} indices: {H_indices}.'
        )

    oh_dipoles = {}
    for i, h_idx in enumerate(H_indices):
        cols = slice(1 + 3 * i, 1 + 3 * i + 3)
        oh_dipoles[h_idx] = data[:, cols]

    return oh_dipoles


def load_OH_types(OH_analysis_file, H_indices):
    """Load the OH type for every selected H from OH_analysis.csv."""
    analysis = pd.read_csv(OH_analysis_file)
    required_columns = {'H_index', 'OH_type'}
    missing_columns = required_columns - set(analysis.columns)
    if missing_columns:
        raise ValueError(
            f'OH analysis file {OH_analysis_file!r} is missing columns '
            f'{sorted(missing_columns)}.'
        )
    if analysis['H_index'].duplicated().any():
        duplicates = analysis.loc[
            analysis['H_index'].duplicated(keep=False), 'H_index'
        ].tolist()
        raise ValueError(
            f'Duplicate H indices in {OH_analysis_file!r}: {duplicates}.'
        )

    OH_types_by_H = dict(zip(analysis['H_index'], analysis['OH_type']))

    missing = [h_idx for h_idx in H_indices if h_idx not in OH_types_by_H]
    if missing:
        raise ValueError(
            f'OH analysis file {OH_analysis_file!r} has no OH type for '
            f'H indices {missing}.'
        )
    return {h_idx: OH_types_by_H[h_idx] for h_idx in H_indices}


def _parse_Hbond_O_indices(value, h_idx, OH_analysis_file):
    """Parse one Hbond_O_indices CSV cell into a list of atom indices."""
    if isinstance(value, str):
        try:
            value = ast.literal_eval(value)
        except (SyntaxError, ValueError) as exc:
            raise ValueError(
                f'Invalid Hbond_O_indices for H index {h_idx} in '
                f'{OH_analysis_file!r}: {value!r}.'
            ) from exc
    if not isinstance(value, (list, tuple)):
        raise ValueError(
            f'Invalid Hbond_O_indices for H index {h_idx} in '
            f'{OH_analysis_file!r}: expected a list, got {value!r}.'
        )
    try:
        return [int(o_idx) for o_idx in value]
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f'Invalid Hbond_O_indices for H index {h_idx} in '
            f'{OH_analysis_file!r}: {value!r}.'
        ) from exc


def load_OH_frequency_metadata(OH_analysis_file, H_indices):
    """Load bond lengths, OH types, and H-bond status for selected Hs."""
    analysis = pd.read_csv(OH_analysis_file)
    required_columns = {
        'H_index', 'OH_type', 'OH_bond_length', 'Hbond_O_indices'
    }
    missing_columns = required_columns - set(analysis.columns)
    if missing_columns:
        raise ValueError(
            f'OH analysis file {OH_analysis_file!r} is missing columns '
            f'{sorted(missing_columns)}.'
        )
    if analysis['H_index'].duplicated().any():
        duplicates = analysis.loc[
            analysis['H_index'].duplicated(keep=False), 'H_index'
        ].tolist()
        raise ValueError(
            f'Duplicate H indices in {OH_analysis_file!r}: {duplicates}.'
        )

    analysis_by_H = analysis.set_index('H_index', drop=False)
    missing = [
        h_idx for h_idx in H_indices
        if h_idx not in analysis_by_H.index
    ]
    if missing:
        raise ValueError(
            f'OH analysis file {OH_analysis_file!r} has no data for '
            f'H indices {missing}.'
        )

    selected = analysis_by_H.loc[H_indices].copy()
    selected['OH_bond_length'] = pd.to_numeric(
        selected['OH_bond_length'], errors='coerce'
    )
    invalid_lengths = selected.loc[
        ~np.isfinite(selected['OH_bond_length']), 'H_index'
    ].tolist()
    if invalid_lengths:
        raise ValueError(
            f'Invalid OH_bond_length in {OH_analysis_file!r} for '
            f'H indices {invalid_lengths}.'
        )

    hbond_indices = [
        _parse_Hbond_O_indices(value, h_idx, OH_analysis_file)
        for h_idx, value in zip(
            selected['H_index'], selected['Hbond_O_indices']
        )
    ]
    selected['Hbond_status'] = [
        'H-bonded' if indices else 'free' for indices in hbond_indices
    ]
    return selected.reset_index(drop=True)


def compute_OH_frequency_statistics(
    freq_cm, ir_by_H, min_freq_cm=OH_STRETCH_MIN_FREQ_CM
):
    """Return intensity-weighted stretching mean and stdev for each H."""
    freq_cm = np.asarray(freq_cm, dtype=float)
    stretch_mask = np.isfinite(freq_cm) & (freq_cm > min_freq_cm)
    if not np.any(stretch_mask):
        raise ValueError(
            f'No finite frequencies above {min_freq_cm:g} cm^-1.'
        )

    stretch_freq = freq_cm[stretch_mask]
    frequency_stats = {}
    for h_idx, intensity in ir_by_H.items():
        intensity = np.asarray(intensity, dtype=float)
        if intensity.shape != freq_cm.shape:
            raise ValueError(
                f'IR intensity for H index {h_idx} has shape '
                f'{intensity.shape}, expected {freq_cm.shape}.'
            )
        weights = intensity[stretch_mask]
        if not np.all(np.isfinite(weights)):
            raise ValueError(
                f'IR intensity for H index {h_idx} contains non-finite '
                f'values above {min_freq_cm:g} cm^-1.'
            )
        if np.any(weights < 0):
            raise ValueError(
                f'IR intensity for H index {h_idx} contains negative '
                f'values above {min_freq_cm:g} cm^-1.'
            )
        total_weight = np.sum(weights)
        if total_weight <= 0:
            raise ValueError(
                f'IR intensity for H index {h_idx} has zero total '
                f'weight above {min_freq_cm:g} cm^-1.'
            )

        mean = np.sum(weights * stretch_freq) / total_weight
        variance = (
            np.sum(weights * (stretch_freq - mean) ** 2) / total_weight
        )
        frequency_stats[h_idx] = (mean, np.sqrt(max(variance, 0.0)))
    return frequency_stats


def write_OH_frequency_data(out_file, frequency_stats):
    """Write weighted stretching-frequency statistics by H index."""
    rows = [
        (h_idx, mean, stdev)
        for h_idx, (mean, stdev) in frequency_stats.items()
    ]
    np.savetxt(
        out_file,
        rows,
        header='H_index OH_freq(cm-1) OH_freq_stdev(cm-1)',
        fmt=['%d', '%.6f', '%.6f'],
    )


def fit_global_frequency_regression(bond_lengths, frequencies):
    """Fit frequency = slope * bond length + intercept and return R^2."""
    bond_lengths = np.asarray(bond_lengths, dtype=float)
    frequencies = np.asarray(frequencies, dtype=float)
    valid = np.isfinite(bond_lengths) & np.isfinite(frequencies)
    bond_lengths = bond_lengths[valid]
    frequencies = frequencies[valid]
    if bond_lengths.size < 2:
        raise ValueError(
            'At least two finite points are required for regression.'
        )
    if np.ptp(bond_lengths) == 0:
        raise ValueError(
            'At least two distinct OH bond lengths are required.'
        )

    slope, intercept = np.polyfit(bond_lengths, frequencies, 1)
    predicted = slope * bond_lengths + intercept
    residual_sum = np.sum((frequencies - predicted) ** 2)
    total_sum = np.sum((frequencies - np.mean(frequencies)) ** 2)
    r_squared = 1.0 if total_sum == 0 else 1.0 - residual_sum / total_sum
    return slope, intercept, r_squared


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


def plot_OH_ir(out_file, freq_cm, ir_by_H, OH_dipole_H_labels=False, freq_range=None, normalize=False):
    plt.figure(figsize=(7, 5))
    for label_idx, (h_idx, ir) in enumerate(ir_by_H.items()):
        if OH_dipole_H_labels is not False:
            label = OH_dipole_H_labels[label_idx]
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


def plot_OH_frequency_vs_bond_length(
    out_file, frequency_stats, frequency_metadata
):
    """Plot OH stretching frequency against bond length by OH category."""
    plot_data = frequency_metadata.copy()
    plot_data['OH_freq'] = [
        frequency_stats[h_idx][0] for h_idx in plot_data['H_index']
    ]
    plot_data['OH_freq_stdev'] = [
        frequency_stats[h_idx][1] for h_idx in plot_data['H_index']
    ]

    fig, ax = plt.subplots(figsize=(7, 5))
    OH_types = list(dict.fromkeys(plot_data['OH_type']))
    colors = plt.get_cmap('tab10').colors
    markers = ('o', 's', '^', 'D', 'v', 'P', 'X', '<', '>')
    for type_idx, OH_type in enumerate(OH_types):
        color = colors[type_idx % len(colors)]
        marker = markers[type_idx % len(markers)]
        for status in ('H-bonded', 'free'):
            group = plot_data[
                (plot_data['OH_type'] == OH_type)
                & (plot_data['Hbond_status'] == status)
            ]
            if group.empty:
                continue
            ax.errorbar(
                group['OH_bond_length'],
                group['OH_freq'],
                yerr=group['OH_freq_stdev'],
                fmt=marker,
                linestyle='none',
                markersize=7,
                markerfacecolor=color if status == 'H-bonded' else 'none',
                markeredgecolor=color,
                ecolor=color,
                capsize=3,
                label=f'{OH_type} ({status})',
            )

    bond_lengths = plot_data['OH_bond_length'].to_numpy()
    frequencies = plot_data['OH_freq'].to_numpy()
    try:
        slope, intercept, r_squared = fit_global_frequency_regression(
            bond_lengths, frequencies
        )
    except ValueError as exc:
        ax.text(
            0.03, 0.97, f'Regression unavailable: {exc}',
            transform=ax.transAxes, va='top'
        )
    else:
        line_x = np.linspace(np.min(bond_lengths), np.max(bond_lengths), 100)
        ax.plot(
            line_x, slope * line_x + intercept,
            color='black', linewidth=1.5, label='Global linear regression'
        )
        ax.text(
            0.03,
            0.97,
            fr'$\nu = {slope:.2f}d_{{OH}} {intercept:+.2f}$' '\n'
            f'$R^2 = {r_squared:.3f}$',
            transform=ax.transAxes,
            va='top',
        )

    ax.set_xlabel(r'O-H bond length ($\mathrm{\AA}$)')
    ax.set_ylabel(r'OH stretching frequency (cm$^{-1}$)')
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_file, dpi=300)
    plt.close(fig)


def get_OH_dipoles_ir(in_file='in.yaml'):
    with open(in_file, 'r') as stream:
        try:
            param = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
            raise

    OH_dipoles_file = param.get('OH_dipoles_out', 'OH_dipoles.dat')
    H_indices = select_H_indices(param, OH_dipoles_file)
    OH_output_plot = param.get('OH_output_plot', 'IR_spectrum_OH.pdf')
    OH_output_group_plot = param.get(
        'OH_output_group_plot', 'IR_spectrum_OH_group.pdf'
    )
    OH_output_group_data = param.get(
        'OH_output_group_data', 'IR_spectrum_OH_group.dat'
    )
    OH_output_data = param.get('OH_output_data', 'IR_spectrum_OH.dat')
    OH_freq_plot = param.get(
        'OH_freq_vs_bond_length_plot', 'OH_freq_vs_bond_length.pdf'
    )
    OH_freq_data = param.get('OH_freq_data', 'OH_freq.dat')
    OH_analysis_file = param.get('OH_analysis_out', 'OH_analysis.csv')
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
    OH_dipole_H_labels = param.get('OH_dipole_H_labels', False)

    oh_dipoles = load_OH_dipoles(OH_dipoles_file, H_indices)

    freq_cm = None
    ir_by_H = {}
    for h_idx in H_indices:
        dipole = oh_dipoles[h_idx]  # (N_frames, 3)
        # [:-1]: keep the length even, matching ir.py's convention
        # (a slightly faster FFT); harmless if the length was already
        # even.
        mu_x, mu_y, mu_z = dipole[:-1].T

        freq_cm_h, ir_h = compute_ir_from_dipole_components(
            [mu_x, mu_y, mu_z], dt_fs,
            use_window=use_window, use_gradient=use_gradient,
            cut_autocorr=cut_autocorr, normalize=normalize,
        )
        if freq_cm is None:
            freq_cm = freq_cm_h
        ir_by_H[h_idx] = ir_h

    write_OH_ir_data(OH_output_data, freq_cm, ir_by_H)
    plot_OH_ir(OH_output_plot, freq_cm, ir_by_H, OH_dipole_H_labels,
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

    if OH_freq_plot is not None or OH_freq_data is not None:
        frequency_stats = compute_OH_frequency_statistics(freq_cm, ir_by_H)
        if OH_freq_data is not None:
            write_OH_frequency_data(OH_freq_data, frequency_stats)
        if OH_freq_plot is not None:
            frequency_metadata = load_OH_frequency_metadata(
                OH_analysis_file, H_indices
            )
            plot_OH_frequency_vs_bond_length(
                OH_freq_plot, frequency_stats, frequency_metadata
            )

    return freq_cm, ir_by_H


if __name__ == '__main__':
    get_OH_dipoles_ir()
