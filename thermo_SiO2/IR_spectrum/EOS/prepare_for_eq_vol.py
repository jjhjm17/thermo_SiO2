#!/usr/bin/env python3
"""This script makes folders for relaxation at the equilibrium volume."""

import os
import sys
import shutil
import numpy as np
from ase.io import read, write
from a_parameters import (
    calculation_folder as calc_folder,
    use_vdw_kernel_file,
    num_samples,
)
from ...util.util import shell
from .util import get_sample_folder_name


def make_input_for_vasp():
    """This function makes input files used by vasp for the EOS calculation."""

    # for debugging
    verbose = True
    # verbose = False

    # print(atoms_and_forces[0]['atoms'].get_positions())
    if os.path.exists("jobList_eq_vol"):
        os.remove("jobList_eq_vol")
    os.chdir(calc_folder)

    eq_vol_relax_folder = "c_eq_vol"
    if os.path.exists(eq_vol_relax_folder):
        print(f"Error: folder {eq_vol_relax_folder} already exists.")
        sys.exit()
    else:
        os.mkdir(eq_vol_relax_folder)
    os.chdir(eq_vol_relax_folder)

    for i_sample in range(num_samples):  # i_sample: index of sample
        sample_folder = get_sample_folder_name(i_sample)
        print(f"\n\n{sample_folder}")
        if not os.path.exists(sample_folder):
            os.mkdir(sample_folder)
        os.chdir(sample_folder)

        root_folder = "../../.."
        template_folder = f"{root_folder}/template"
        snapshot = read(
            f"{root_folder}/{calc_folder}/a_given_vol_relax/{sample_folder}/CONTCAR"
        )

        eq_vol = np.loadtxt(
            f"{root_folder}/{calc_folder}_result/{sample_folder}/EVinet"
        )[
            1
        ]  # Ang^3/atom
        num_atoms = len(snapshot)
        alat = (eq_vol * num_atoms) ** (1 / 3)
        print(f"{alat = :.3f} Ang")
        snapshot.set_cell([[alat, 0, 0], [0, alat, 0], [0, 0, alat]], scale_atoms=True)
        if verbose:
            print("snapshot = ", snapshot)

        write("POSCAR", snapshot, direct=True, vasp5=True)
        shutil.copy(f"{template_folder}/INCAR", ".")
        shutil.copy(f"{template_folder}/INCAR.preconverge.change", ".")
        shutil.copy(f"{template_folder}/KPOINTS", ".")
        os.symlink(f"{template_folder}/POTCAR", "POTCAR")
        if use_vdw_kernel_file:
            os.symlink(f"{template_folder}/vdw_kernel.bindat", "vdw_kernel.bindat")
        shell(f"pwd >> {root_folder}/jobList_eq_vol")
        os.chdir("..")

    os.chdir("..")


if __name__ == "__main__":
    make_input_for_vasp()
