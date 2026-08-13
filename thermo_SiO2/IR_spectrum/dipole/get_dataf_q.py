"""
Convert a POSCAR file to the
'charge' atom_style, using ASE.

Reads:  POSCAR
Writes: coords_q.dataf (atom_style = charge)

Notes
-----
- 'atomic' style datafiles don't carry per-atom charge, so all atoms are
  written out with charge = 0.0. Edit the `charges` dict below (or the
  written file) if you need non-zero starting charges per element.
- The atom-type -> element mapping (Si, O, H, Al = types 1, 2, 3, 4) must
  match the order used when the original data file was created.
"""

from ase.data import atomic_numbers
from ase.io.lammpsdata import read_lammps_data, write_lammps_data
from ase.io import read
from thermo_SiO2.util.util import read_in_yaml

def get_dataf_q_from_POSCAR():

    param = read_in_yaml()
    INPUT_FILE = param['born_poscar']
    OUTPUT_FILE = param['coord_file']

    # Atom type (as used in the LAMMPS data file) -> chemical species, in order.
    SPECORDER = param['atom_symbols'].split()
    Z_OF_TYPE = {i + 1: atomic_numbers[sym] for i, sym in enumerate(SPECORDER)}

    # Optional: per-element starting charges (all zero by default).
    # CHARGES = {"Si": 0.0, "O": 0.0, "H": 0.0, "Al": 0.0}
    CHARGES = {key:0.0 for key in SPECORDER}


    # atoms = read_lammps_data(
    #     INPUT_FILE,
    #     Z_of_type=Z_OF_TYPE,
    #     atom_style="atomic",
    #     sort_by_id=True,
    # )
    atoms = read(INPUT_FILE)


    charges = [CHARGES[sym] for sym in atoms.get_chemical_symbols()]
    atoms.set_initial_charges(charges)

    write_lammps_data(
        OUTPUT_FILE,
        atoms,
        specorder=SPECORDER,
        atom_style="charge",
        masses=True,
    )

    print(f"Read {len(atoms)} atoms from '{INPUT_FILE}'")
    print(f"Wrote '{OUTPUT_FILE}' with atom_style='charge'")


if __name__ == "__main__":
    get_dataf_q_from_POSCAR()

