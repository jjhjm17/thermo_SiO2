"""This function run_neb collects some tools for neb.
1) lammps like output (reaction coordinates, their total length)
2) regular + CI neb
3) Repel initial interpolated images using the Lennard-Jones potential.
4) Interpolate from previous xyz file even when the number of images is different in a piecewise-linear way.


This script performs the nudged elastic band calculation (NEB) of
(for example) 8 images, with the initial (for example) ../a.a.initial/cfg_d_relaxed.dataf
and the final (for example) ../b.a.final/cfg_d_relaxed.dataf.
(The number 8 includes endpoints.)
The convergence criterion is (for example) 0.01 eV/Ang.
We use grace_model (for example) 
'/home/atuin/a102cb/a102cb12/SiO2_mesoporous/b.MTP/n.grace/a.2L_medium/a.medium_5Ang/seed/1/model.yaml'

0) First use repelsive LJ potential to repeal atoms whose bond lengths are < repel_cutoff when first_repel_LJ = True

1) We do NEB

2) Like lammps, do climbing NEB.

Use all images in 1 GPU.

If lmp_neb_dump is present, then use them as initial images.
Since the initial and final images are relaxed using GRACE, when we read dump files in lmp_neb_dump, read dump.2 (1-based), dump.3, ... dump.7 (the second last one) and use them as guesses for the 2nd to 2nd last images.
The default behavior is not to read guess.xyz (see use_interpolate_xyz, interpolate_xyz variable).
When interpolate_xyz_file's number of images is different from n_images, do the piecewise linear interpolation.

Ref https://ase.gitlab.io/ase/examples_generated/03-tutorials/neb_selfdiffusion.html

During neb minimization, could we print out at the screen information like lammps neb. For now we need only Step, MaxReplicaForce (fmax in min.log), RDT, RD1, PE1, .., RD8, PE8 (the last column.) with interval step_output 

   Setting up regular NEB ...
       Step     MaxReplicaForce  MaxAtomForce      GradV0         GradV1         GradVc          EBF            EBR            RDT            RD1            PE1            RD2            PE2            RD3            PE3            RD4            PE4            RD5            PE5            RD6            PE6            RD7            PE7            RD8            PE8       
            0   32.782269        298.17659      0.0099121631   0.0083163153   32.795154      3.6671868      4.323187       3.5660299      0              -7299.4247     0.14285714     -7298.5941     0.28571429     -7296.8917     0.42857143     -7295.7575     0.57142857     -7296.0675     0.71428571     -7297.5833     0.85714286     -7299.2852     1              -7300.0807     
   
   ...
   Setting up climbing ...
   Climbing replica = 1
       Step     MaxReplicaForce  MaxAtomForce      GradV0         GradV1         GradVc          EBF            EBR            RDT            RD1            PE1            RD2            PE2            RD3            PE3            RD4            PE4            RD5            PE5            RD6            PE6            RD7            PE7            RD8            PE8       
        13357   0.0099987123     1.3714391e-06  1.1761583e-05  2.1831487e-06  1.1761583e-05  0              0.8171387      15.03718       0              -7299.4255     0.14195051     -7299.6068     0.28517594     -7299.6274     0.42813387     -7299.6318     0.57001136     -7299.6318     0.71406298     -7300.0575     0.85706515     -7300.2424     1              -7300.2426     
   
Keep this documentation as it is hand-written.
"""

from contextlib import ExitStack
from numbers import Integral
import os.path

import numpy as np
from ase.optimize import FIRE
from ase.mep import NEB, NEBTools
from ase.geometry import find_mic
from ase.utils.forcecurve import fit_images
from ase.units import fs
from ase.io import write, read, Trajectory
from ase.calculators.lj import LennardJones

def resample_images(images, n_images):
    """Linearly resample a path over image indices, accounting for periodicity."""
    if n_images < 2:
        raise ValueError('n_images must include at least two endpoints.')
    if len(images) < 2:
        raise ValueError('The interpolation XYZ must contain at least two images.')

    reference = images[0]
    for image_number, image in enumerate(images[1:], start=2):
        if len(image) != len(reference):
            raise ValueError(
                f'XYZ image {image_number} contains {len(image)} atoms; '
                f'expected {len(reference)}.'
            )
        if not np.array_equal(image.numbers, reference.numbers):
            raise ValueError(
                f'Atom ordering in XYZ image {image_number} does not match '
                'the first image.'
            )
        if not np.array_equal(image.pbc, reference.pbc):
            raise ValueError(
                f'Periodicity in XYZ image {image_number} does not match '
                'the first image.'
            )
        if not np.allclose(image.cell, reference.cell, rtol=0, atol=1e-8):
            raise ValueError(
                f'Cell in XYZ image {image_number} does not match the first image.'
            )

    if len(images) == n_images:
        return images

    # Unwrap adjacent source images before sampling across multiple segments.
    positions = [reference.positions.copy()]
    for previous, current in zip(images[:-1], images[1:]):
        displacement, _ = find_mic(
            current.positions - previous.positions, reference.cell, reference.pbc
        )
        positions.append(positions[-1] + displacement)

    resampled = [reference.copy()]
    for coordinate in np.linspace(0, len(images) - 1, n_images)[1:-1]:
        segment = int(coordinate)
        fraction = coordinate - segment
        image = images[segment].copy()
        image.set_positions(
            positions[segment]
            + fraction * (positions[segment + 1] - positions[segment]),
            apply_constraint=False,
        )
        resampled.append(image)
    resampled.append(images[-1].copy())
    return resampled


def print_neb_status(neb_object, optimizer, step_offset=0):
    """Print ASE fmax, path length, and each image's coordinate and energy."""
    band_images = neb_object.images
    if optimizer.nsteps == 0:
        headers = ['MaxReplicaForce', 'RDT']
        for image_number in range(1, len(band_images) + 1):
            headers.extend([f'RD{image_number}', f'PE{image_number}'])
        print(f'{"Step":>10}' + ''.join(f'{name:>16}' for name in headers),
              flush=True)

    # Match ASE's per-atom fmax in min.log, rather than LAMMPS's replica norm.
    forces = neb_object.get_forces()
    fmax = np.linalg.norm(forces, axis=1).max()

    # Cumulative 3N-dimensional distances, using the same MIC as ASE NEB.
    rd = [0.0]
    for previous, current in zip(band_images[:-1], band_images[1:]):
        displacement, _ = find_mic(
            current.positions - previous.positions, previous.cell, previous.pbc
        )
        rd.append(rd[-1] + np.linalg.norm(displacement))
    rdt = rd[-1]
    rd = np.array(rd)
    if rdt > 0:
        rd /= rdt

    values = [fmax, rdt]
    for coordinate, image in zip(rd, band_images):
        values.extend([coordinate, image.get_potential_energy()])
    step = step_offset + optimizer.nsteps
    print(f'{step:10d}' + ''.join(f'{value:16.8g}' for value in values),
          flush=True)


def get_FIRE(neb_object):
    # return FIRE(neb_object, logfile='min.log', dt=0.5*fs, dtmax=2*fs)
    # return FIRE(neb_object, logfile='min.log', dt=0.25*fs, dtmax=1*fs)
    return FIRE(neb_object, logfile='min.log', dt=0.1*fs, dtmax=1*fs, maxstep=0.02)
     # see  https://henkelmangroup.github.io/vtsttools/optimizers.html
    # dt: initial time step, 0.5 fs: for hydrogen
    # dt max: 10 times MD step
    #   "Thus the only adjustable parameter of the method is the maximum time step tmax."
    #   FIRE paper, https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.97.170201

def run_neb(
    cfg_initial,
    cfg_final,
    grace_dir,
    *,
    n_images=8,
    k_spring=0.1,
    ftol_normal_neb=0.05,
    ftol_climb_neb=0.05,
    symbols='Si O H Al',
    lmp_neb_dump=None,
    interpolate_xyz_file=None,
    dump_interval=10,
    step_output=10,
    n_max_steps=10000,
    repel_cutoff=0.8,
    debug_lj=False,
):
    """Run regular and climbing-image NEB and return the final ASE images.

    ``cfg_initial`` and ``cfg_final`` are ASE Atoms objects. ``n_images``
    includes both endpoints; ``n_max_steps`` applies to each NEB stage.
    ``grace_dir`` is the TensorPotential model directory and is unused when
    ``debug_lj=True``. ``symbols`` maps LAMMPS atom types when reading dumps.

    Initial paths are selected in this order: LAMMPS intermediate dumps,
    an interpolation XYZ when ``interpolate_xyz_file`` is not None, then
    linear interpolation between the supplied endpoints. An XYZ supplies
    the entire path, including its endpoints, and is resampled when its
    image count differs. A specified XYZ must exist; use None to generate
    a linear path.

    All relative input paths and output filenames use the caller's current
    working directory. Existing outputs are overwritten, except trajectories
    are appended during climbing NEB after being created during regular NEB.
    Input Atoms objects may be modified and assigned calculators.

    ``k_spring`` is passed as ``k`` to both NEB stages; its default is
    0.1 eV/Ang^2. LJ repulsion always runs on intermediate images, preserving
    the original behavior of the unused first_repel variable.
    """
    for name, value, minimum in (
        ('n_images', n_images, 3),
        ('step_output', step_output, 1),
        ('dump_interval', dump_interval, 1),
        ('n_max_steps', n_max_steps, 0),
    ):
        if (isinstance(value, bool) or not isinstance(value, Integral)
                or value < minimum):
            raise ValueError(f'{name} must be an integer >= {minimum}.')

    # Keep imports free of model loading and optional TensorFlow dependencies.
    import matplotlib.pyplot as plt
    from thermo_SiO2.io import read_sil
    if not debug_lj:
        from tensorpotential.calculator import TPCalculator

    with ExitStack() as outputs:
        if lmp_neb_dump:
            dump_by_image = {
                int(os.path.basename(dump).rsplit('.', 1)[1]): dump
                for dump in lmp_neb_dump
            }
            intermediate_image_numbers = range(2, n_images)
            missing_dumps = [
                image_number
                for image_number in intermediate_image_numbers
                if image_number not in dump_by_image
            ]
            if missing_dumps:
                raise FileNotFoundError(
                    f'Missing LAMMPS NEB dump files for images: {missing_dumps}'
                )

            print('\nLAMMPS NEB dumps exist. We use their final snapshots as intermediate images.\n')
            images = [cfg_initial]
            for image_number in intermediate_image_numbers:
                dump_images = read_sil(
                    dump_by_image[image_number], index=-1, atom_symbols=symbols
                )
                image = dump_images[0]
                if len(image) != len(cfg_initial):
                    raise ValueError(
                        f'{dump_by_image[image_number]} contains {len(image)} atoms; '
                        f'expected {len(cfg_initial)}.'
                    )
                if image.get_chemical_symbols() != cfg_initial.get_chemical_symbols():
                    raise ValueError(
                        f'Atom ordering in {dump_by_image[image_number]} does not match '
                        'the GRACE-relaxed endpoint.'
                    )
                images.append(image)
            images.append(cfg_final)
            interpolate_images = False
        elif interpolate_xyz_file is not None:
            interpolate_images = False
            print(f'\nPrevious interpolation exists. We read it: {interpolate_xyz_file}.\n')
            # It takes time (some minutes) to guess interpolation.
            images = read(interpolate_xyz_file, index=':')
            source_count = len(images)
            images = resample_images(images, n_images)
            if source_count != n_images:
                print(
                    f'Linearly resampled {interpolate_xyz_file}: '
                    f'{source_count} -> {n_images} images.\n'
                )
        else:
            interpolate_images = True
            print('\nPrevious interpolation does not exists. We make one.\n')
            images = [cfg_initial]
            for i in range(n_images - 2):
                images.append(cfg_initial.copy())
            images.append(cfg_final)

        if len(images) != n_images:
            raise ValueError(f'Expected {n_images} NEB images, but loaded {len(images)}.')



        neb = NEB(images, k=k_spring, method='improvedtangent',
                  allow_shared_calculator=True)
        # improvedtangent: default
        if interpolate_images:
            # neb.interpolate(mic=True, method='idpp')
            # NEB with IDPP: Image Dependent Pair Potential for improved interpolation of NEB initial guess
            # Do not use IDPP, H is too far away.
            neb.interpolate(mic=True)  # linear guess

            write('interpolate.xyz', images, format='extxyz')



        # repel atoms < repel_cutoff
        sigma = repel_cutoff / 2**(1/6)

        LJrepulsive = LennardJones(
            epsilon=0.1,   # eV
            sigma=sigma,
            rc=repel_cutoff,
        )
        for i_image in range(1, n_images-1):
            print(f'\n\n{i_image = }')
            image = images[i_image]
            image.calc = LJrepulsive
            dyn = FIRE(image, maxstep=0.01)  # maxstep: Ang
            dyn.run(fmax=0.01)

        write('interpolate_repel.xyz', images, format='extxyz')

        # if debug_lj:
        #     for image in images:
        #         image.calc = LennardJones(sigma=1.5, epsilon=1)
        # else:
        #     calc = TPCalculator(grace_dir,
        #                         # pad_neighbors_fraction = 0.05,
        #                         pad_neighbors_fraction = 0.1,
        #                         pad_atoms_number = 2,
        #                         # min_dist=0.5) # initial guess 0.45
        #                         min_dist=0.4)
        #     use_shared_calc = False
        #     if use_shared_calc:
        #         neb.set_calculators(calc)
        #     else:
        #         neb.set_calculators([calc] * n_images)

        for image in images:
            # breakpoint()
            # image.calc = calc
            # separate calculator objects to reuse previous calculations
            # test case: shared calc: 7 sec / step, not shared: 2 sec/step
            # system: 700 atoms, 12 images

            if debug_lj:
                image.calc = LennardJones(sigma=1.5, epsilon=1)
            else:
                image.calc = TPCalculator(grace_dir,
                                 # pad_neighbors_fraction = 0.05,
                                 pad_neighbors_fraction = 0.1,
                                 pad_atoms_number = 2,
                                 # min_dist=0.5) # initial guess 0.45
                                 min_dist=0.4)


        pot_E_step_0= [image.get_potential_energy() for image in images]
        print(f'{pot_E_step_0 = }')

        # minimizer = MDMin(neb, logfile='min.log', maxstep=0.01)  # oscillate fmax 0.7 0.5 0.7 ..
        # maxstep: Ang
        # minimizer = FIRE(neb, logfile='min.log', maxstep=0.05, dt=0.05)
        # dt: initial time step
        # minimizer = FIRE(neb, logfile='min.log')  # diverges

        print('\nSetting up regular NEB ...', flush=True)
        minimizer = get_FIRE(neb)
        minimizer.attach(print_neb_status, interval=step_output,
                         neb_object=neb, optimizer=minimizer)

        for i in range(1, n_images-1):
            minimizer.attach(outputs.enter_context(
                Trajectory('neb-%d.traj' % i, 'w', images[i])),
                            interval=dump_interval)
        minimizer.run(fmax=ftol_normal_neb, steps=n_max_steps)
        if minimizer.nsteps % step_output:
            print_neb_status(neb, minimizer)
        regular_neb_steps = minimizer.nsteps

        # see https://docs.ase-lib.org/examples_generated/03-tutorials/diffusion_neb.html#diffusion-tutorial
        fit = fit_images(images)
        path = np.array(fit.path)
        energies = np.array(fit.energies)
        path /= path[-1]
        print(f'{path =}')      # reaction coordinate for each NEB image
        print(f'{energies =}')  # relative energies

        write('neb_a_regular.xyz', images, format='extxyz')
        outputs.close()

        # climbing image
        print('\nSetting up climbing NEB ...', flush=True)
        nebClimb = NEB(images, k=k_spring, method='improvedtangent',
                       allow_shared_calculator=True, climb=True)
        # minimizer = FIRE(nebClimb, logfile='min.log', maxstep=0.05, dt=0.05)
        # minimizer = FIRE(nebClimb, logfile='min.log')
        minimizer = get_FIRE(nebClimb)
        minimizer.attach(print_neb_status, interval=step_output,
                         neb_object=nebClimb, optimizer=minimizer,
                         step_offset=regular_neb_steps)
        for i in range(1, n_images-1):
            minimizer.attach(outputs.enter_context(
                Trajectory('neb-%d.traj' % i, 'a', images[i])),  # a: append
                            interval=dump_interval)
        minimizer.run(fmax=ftol_climb_neb, steps=n_max_steps)
        if minimizer.nsteps % step_output:
            print_neb_status(nebClimb, minimizer, step_offset=regular_neb_steps)

        nebtools = NEBTools(images)
        E_barrier, delta_E = nebtools.get_barrier(fit=False)  # fit=False: no interpolation between highest images

        # see https://docs.ase-lib.org/examples_generated/03-tutorials/diffusion_neb.html#diffusion-tutorial
        fit = fit_images(images)
        path = np.array(fit.path)
        energies = np.array(fit.energies)
        path /= path[-1]
        print(f'{path =}')      # reaction coordinate for each NEB image
        print(f'{energies =}')  # relative energies


        fig, ax = plt.subplots()
        ax.plot(
            path, energies, marker='+',
        )
        ax.set_xlabel('Reaction coordinate (normalized)')
        ax.set_ylabel('Energy (eV)')
        fig.savefig('fig_neb.pdf', transparent=True)
        plt.close(fig)

        print(f'The energy barrier is {E_barrier:.4f} eV.')
        print(f'delta E: {delta_E:.4f} eV.')

        write('neb_b_climb.xyz', images, format='extxyz')
        return images
