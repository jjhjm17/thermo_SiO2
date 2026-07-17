"""This script obtains the diffusion coefficient of O atoms from vasp OUTCAR."""
import numpy as np
from ase.io import read
from ase.md.analysis import DiffusionCoefficient
from ase.units import fs
import a_parameters as param
import matplotlib.pyplot as plt


def unwrap(images):
    """This function unwraps the coordinates in images, or ase trajectories."""
    cell = images[0].cell
    is_diagonal = np.allclose(cell, np.diag(np.diagonal(cell)), atol=1e-6)
    if not is_diagonal:
        raise ValueError('Since the cell is not diagonal, the unwrapping is not yet implemented.')
    positions = [atoms.get_positions() for atoms in images]
    unwrapped = [positions[0].copy()]

    for i in range(1, len(positions)):
        delta = positions[i] - positions[i-1]
        # detect crossings (assuming orthorhombic)
        delta -= np.round(delta / cell.diagonal()) * cell.diagonal()
        # test
        # delta =  0.1, diagonal =  1, delta -= np.round(0.1/1)*1  -> -=  0
        # delta =  1.1, diagonal =  1, delta -= np.round(1.1/1)*1  -> -=  1
        # delta = -0.1, diagonal =  1, delta -= np.round(-0.1/1)*1 -> -=  0
        # delta = -1.1, diagonal =  1, delta -= np.round(-1.1/1)*1 -> -= -1, += 1
        unwrapped.append(unwrapped[-1] + delta)
        images[i].positions = unwrapped[i]
    return images


def write_center_of_mass(traj, filename):
    """Write trajectory centered at center of mass to a file."""
    with open(filename, 'w') as f:
        f.write(f"# Frame   positions of the center of mass (Ang)\n")
        for i_traj, atoms in enumerate(traj):
            com = atoms.get_center_of_mass()
            # centered_positions = atoms.get_positions() - com
            f.write(f"{i_traj}  {com[0]:.6f} {com[1]:.6f} {com[2]:.6f}\n")
    print(f"Center of mass trajectory written to {filename}")


def get_diffusion_coeff_traj(traj, step, show_msd_plot=True, run_idx=None):
    """Get diffusion coefficient of O atoms from ase trajectory of vasp."""
    diff_timestep = step * fs

    # if hasattr(param, 'unwrap') and param.unwrap:
    traj = unwrap(traj)
    # vasp trajectory is always wrapped, so unwrapping is needed.
    if hasattr(param, 'write_center_of_mass') and param.write_center_of_mass:
        # write centered trajectory
        write_center_of_mass(traj, f'center_of_mass_{run_idx}.txt')

    # --- Select only O atoms ---
    # Assuming traj is a list of Atoms objects
    atoms0 = traj[0]
    O_indices = [atom.index for atom in atoms0 if atom.symbol == 'O']

    # --- Slice trajectory from step_start onwards and subsample every n_every ---
    # traj_sub = traj[step_start::n_every]

    # --- Compute diffusion coefficient ---
    diff = DiffusionCoefficient(traj, timestep=diff_timestep, atom_indices=O_indices)

    # Perform the actual calculation
    slopes, std = diff.get_diffusion_coefficients()
    slope = slopes[0] * (1000 * fs)  # Ang^2/(ase time unit) -> Ang^2 / ps
    std = std[0] * (1000 * fs)  # unit  Ang^2/ps
    # breakpoint()
    print(f"Diffusion coefficient (O atoms): {slope:.4e} +- {std:.4e} (std) Ang^2/ps")

    # Optional: get MSD and plot
    # msd, times = diff.get_msd()
    # plt.plot(times, msd)
    # plt.xlabel('Time (fs)')
    # plt.ylabel('MSD (Å²)')
    # plt.title('Mean Square Displacement (O atoms)')
    # plt.savefig('fig__diffuse.pdf')
    # plt.show()
    diff.print_data()
    # diff.plot(show=True)
    diff.plot(show=show_msd_plot)
    return slope, std


def get_diffusion_coeff():
    # --- Parameters ---
    file = param.file          # path to your OUTCAR
    # interval = param.interval
    index = param.index
    # --- Read trajectory ---
    # ASE can read multiple frames from OUTCAR
    traj = read(file, index=index)  # read all snapshots
    step = param.step

    slope, std = get_diffusion_coeff_traj(traj, step)
    return slope, std

if __name__ == "__main__":
    get_diffusion_coeff()
