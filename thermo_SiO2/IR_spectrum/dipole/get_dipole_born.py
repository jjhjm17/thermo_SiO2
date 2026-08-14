"""This function obtains the dipole moment from the Born effective charge tensor and lammps displacement.

Input
in.yaml

output
dipole_out
"""
import numpy as np
import yaml
from ase.io import read
from ase.geometry import find_mic
from thermo_SiO2.io import read_sil


def read_born_charges(born_file, cfg_0=None):
    """Parse a BORN file: line 1 is a comment, line 2 is epsilon (ignored),
    remaining lines are one 3x3 Born effective charge tensor per atom
    (9 values: xx, xy, xz, yx, yy, yz, zx, zy, zz).
    https://phonopy.github.io/phonopy/input-files.html#born-file
    """
    with open(born_file, 'r') as f:
        lines = f.readlines()

    # line 0: comment ("# epsilon and Z* of atoms ...")
    # line 1: epsilon tensor -> ignored
    data_lines = lines[2:]

    Z = np.array([[float(x) for x in line.split()] for line in data_lines])
    n_atoms = Z.shape[0]
    Z = Z.reshape(n_atoms, 3, 3)  # shape (N_atoms, 3, 3)
    # DEBUG = True
    # if DEBUG:
    #     print(f'{Z[:3] = }')

    # DEBUG_PARTIAL = True
    DEBUG_PARTIAL = False
    if DEBUG_PARTIAL:
        print('DEBUG_PARTIAL')
        # DEBUG_H_ONLY = True
        DEBUG_H_ONLY = False
        if DEBUG_H_ONLY:
            partial_atoms = ['H']
            print('DEBUG_H_ONLY')
        else:
            # partial_atoms = ['O']
            # partial_atoms = ['Si']
            partial_atoms = ['Al']
            print(f'{partial_atoms = }')
        for i_atom, symbol in enumerate(cfg_0.symbols):
            if symbol not in partial_atoms:
                Z[i_atom, :, : ] = np.zeros(3)

    return Z

def get_fixed_Z(cfg_0):
    nominal_charges = {'Si': 1.2, 'O': -0.6, 'H': 0.3, 'Al': 0.9}

    atoms = cfg_0
    n_atoms = len(cfg_0)
    Z = np.array([np.eye(3) * nominal_charges[s] for s in atoms.symbols])
    return Z


def write_dipoles(dipole_out, dipoles):
    with open(dipole_out, 'w') as f:
        f.write('# TimeStep dipole moment  x, y, z (|e|)\n')
        for i, dipole in enumerate(dipoles):
            f.write('{} {:.5f} {:.5f} {:.5f}\n'.format(i, *dipole))


def get_dipole_born(in_file='in.yaml'):
    with open(in_file, 'r') as stream:
        try:
            param = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    cfgs = read_sil(param['dump_unfolded'],
                    atom_symbols=param['atom_symbols'])
    print('cfgs were read.')

    
    if not param.get('test_fixed_charge_Si_O_H_Al', False):
        # Check of born_poscar and dump file have the same atomic 
        # structures, and the ordering of atomic symbols are the same.
        cfg_born = read(param['born_poscar'])
        cfg0 = cfgs[0]
        if not (cfg0.symbols == cfg_born.symbols).all():
            # .all(): if all values are true
            raise ValueError('The atomic symbols of born_poscar and the dump file are differernt.')
        elif not (cfg0.cell == cfg_born.cell).all():
             # .all(): if all values are true
             raise ValueError('The cells of born_poscar and the dump file are differernt.')

        dr = cfg0.positions - cfg_born.positions 
        _, dr_dist = find_mic(dr, cell=cfg0.cell, pbc=cfg0.pbc)
        dr_max = np.max(dr_dist)
        print(f'The max displacement between the positions of born_poscar and dump file: {dr_max:.2f} Ang.')

        if 5 > dr_max > 3:
            print(f'\nWarning: The max displacement between the positions of born_poscar and dump file is large, {dr_max:.2f} Ang. please check if the structures are the same, and the atomic order is not changed.\n') 
        elif dr_max >= 5:
            raise ValueError(f'The max displacement between the positions of born_poscar and dump file is very large, {dr_max:.2f} Ang. Please check if the structures are the same, and the atomic order is not changed.') 


    dipole_out = param.get('dipole_out')

    # if param.get('test_fixed_charge_Si_O_H_Al', False):
    #     # https://doi.org/10.1063/5.0194486
    #     nominal_charges = {'Si': 1.2, 'O': -0.6, 'H': 0.3, 'Al': 0.9}

    #     dipoles = []
    #     for atoms in cfgs:
    #         symbols = atoms.get_chemical_symbols()
    #         positions = atoms.get_positions()
    #         charges = np.array([nominal_charges[s] for s in symbols])  # (N,)
    #         dipole = (charges[:, None] * positions).sum(axis=0)        # (3,)
    #         dipoles.append(dipole)

    #     if dipole_out is not None:
    #         write_dipoles(dipole_out, dipoles)
    #     return dipoles

    
    if param.get('test_fixed_charge_Si_O_H_Al', False):
        Z = get_fixed_Z(cfg_0 = cfgs[0])
    else:
        Z = read_born_charges(param['born_file'], cfg_0 = cfgs[0])  # shape (N_atoms, 3, 3)

    print(f'total {len(cfgs)} cfgs')
    dipoles = []
    for i, atoms in enumerate(cfgs):
        if i % 1000 == 0:
            print(f'{i} / {len(cfgs)} cfgs')
        positions = atoms.get_positions()  # shape (N_atoms, 3)

        if positions.shape[0] != Z.shape[0]:
            raise ValueError(
                f"Number of atoms in config ({positions.shape[0]}) "
                f"does not match number of Born tensors ({Z.shape[0]})"
            )

        # vasp Z*_k,ij = Omega / e round P_i / round u_k,j(q=0)
        # https://vasp.at/wiki/Born_effective_charges
        #
        # sum_n (Z_n @ r_n) -> total dipole vector, shape (3,)
        dipole = np.einsum('nij,nj->i', Z, positions)
        dipoles.append(dipole)

    if dipole_out is not None:
        write_dipoles(dipole_out, dipoles)
    return dipoles


if __name__ == '__main__':
    get_dipole_born()
