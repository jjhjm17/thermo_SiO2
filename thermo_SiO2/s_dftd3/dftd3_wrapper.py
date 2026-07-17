import subprocess
import os
import json
from numpy import array
from ase.calculators.calculator import Calculator, all_changes
from ase import units

class SDftd3Calculator(Calculator):
    """We use a compiled fortran binary s-dftd3, which is faster than the provided
    python interface."""
    implemented_properties = ["energy", "forces"]

    def __init__(self, workdir=None, **kwargs):
        super().__init__(**kwargs)

        # Use $SCRATCH if available for fast I/O
        if workdir is None:
            # shm_path = "/dev/shm"
            shm_path = os.environ['SCRATCH']
            if os.path.isdir(shm_path) and os.access(shm_path, os.W_OK):
                workdir = os.path.join(shm_path, "dftd3_tmp")
            else:
                workdir = "dftd3_tmp"  # fallback

        self.workdir = workdir
        os.makedirs(self.workdir, exist_ok=True)

    def calculate(self, atoms=None, properties=["energy"], system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)

        # periodic boundary condition
        if atoms.get_pbc().all():
            file_name = 'POSCAR'
        else:
            file_name = 'mol.xyz'
        structure_file = os.path.join(self.workdir, file_name)
        atoms.write(structure_file)

        # call s-dftd3
        subprocess.run(
            ["s-dftd3", file_name, "--zero", "pbe", "--grad", "--json", "--noedisp"],
            cwd=self.workdir,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            text=True
        )

        # parse output (you need to implement this)
        with open(f"{self.workdir}/dftd3.json") as f:
            json_out = json.load(f)
        energy = json_out["energy"]
        forces = - array(json_out['gradient'], dtype=float).reshape(-1, 3)

        # remove files
        os.remove(structure_file)
        os.remove(f"{self.workdir}/dftd3.json")

        # unit conversion from the atomic units to eV, eV/Ang
        energy *= units.Hartree
        forces *= units.Hartree / units.Bohr

        self.results = {
            "energy": energy,
            "forces": forces
        }

