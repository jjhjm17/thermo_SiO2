#!/usr/bin/env python3
"""Generate deterministic synthetic O-H vectors for the IR unit tests.

Each vector oscillates in one normal mode around a fixed equilibrium
vector. The wavenumbers and timestep are chosen to match fixtures/in.yaml.
"""
import argparse
from pathlib import Path

import numpy as np


LIGHT_SPEED_CM_PER_FS = 2.99793e-5
DT_FS = 0.5
N_FRAMES = 4000

# (H index, wavenumber in cm^-1, equilibrium vector in Ang,
#  oscillation-amplitude vector in Ang)
OH_MODES = (
    (700, 3600.0, (0.90, 0.20, 0.10), (0.080, 0.024, 0.016)),
    (682, 3800.0, (-0.85, 0.25, -0.15), (0.040, 0.008, 0.012)),
)


def generate_OH_dipoles():
    """Return steps and synthetic vectors with shape (frames, Hs, xyz)."""
    steps = np.arange(N_FRAMES)
    time_fs = steps * DT_FS
    vectors = []
    for _, wavenumber, equilibrium, amplitude in OH_MODES:
        phase = (
            2.0 * np.pi * wavenumber * LIGHT_SPEED_CM_PER_FS * time_fs
        )
        vector = (
            np.asarray(equilibrium)[None, :]
            + np.cos(phase)[:, None] * np.asarray(amplitude)[None, :]
        )
        vectors.append(vector)
    return steps, np.stack(vectors, axis=1)


def write_OH_dipoles(output_file):
    """Write OH_dipoles.dat in the format consumed by get_OH_dipoles_ir."""
    steps, vectors = generate_OH_dipoles()
    H_indices = [mode[0] for mode in OH_MODES]
    header_columns = ' '.join(
        f'H_{h_idx}_x H_{h_idx}_y H_{h_idx}_z' for h_idx in H_indices
    )
    columns = np.column_stack((steps, vectors.reshape(N_FRAMES, -1)))
    np.savetxt(
        output_file,
        columns,
        header=(
            f'step {header_columns}  (OH dipole components, |e| Ang)'
        ),
        fmt=['%d'] + ['%.5f'] * (3 * len(H_indices)),
    )


def main():
    default_output = Path(__file__).with_name('OH_dipoles.dat')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--output', type=Path, default=default_output,
        help=f'output path (default: {default_output})',
    )
    args = parser.parse_args()
    write_OH_dipoles(args.output)
    frequencies = ', '.join(
        f'H {h_idx}: {wavenumber:g} cm^-1'
        for h_idx, wavenumber, _, _ in OH_MODES
    )
    print(f'Wrote {N_FRAMES} frames to {args.output} ({frequencies})')


if __name__ == '__main__':
    main()
