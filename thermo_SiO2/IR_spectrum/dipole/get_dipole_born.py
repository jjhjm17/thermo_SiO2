"""This function obtains the dipole moment from the Born effective charge tensor and lammps displacement.

Input
in.yaml
        # see explanation of get_OH_dipoles.py for more details
        dump_unfolded : '../a.traj/config.dump'
        atom_symbols : 'Si O H Al'
        charge : 'formal'  # or 'born_isotropic' or 'born_full'
        born_file : 'xxx/sample_0_BORN'  # needed for 'born_isotropic' or 'born_full'
        born_poscar : 'xxx/POSCAR'       # needed for 'born_isotropic' or 'born_full'

output
        born_dipole_out : 'born_dipole.out'
"""
import numpy as np
import yaml
from thermo_SiO2.IR_spectrum.dipole.config_paths import resolve_config_path
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
    formal_charges = {'Si': 1.2, 'O': -0.6, 'H': 0.3, 'Al': 0.9}

    atoms = cfg_0
    n_atoms = len(cfg_0)
    Z = np.array([np.eye(3) * formal_charges[s] for s in atoms.symbols])
    return Z


def get_charge_tensors(param, cfg0):
    """Return charge tensors using the same charge modes as OH analysis."""
    if param.get('charge') == 'formal':
        return get_fixed_Z(cfg0)
    from thermo_SiO2.IR_spectrum.dipole.get_OH_dipoles import (
        get_charge_tensors as get_oh_charge_tensors,
    )
    return get_oh_charge_tensors(param, cfg0)


def write_dipoles(born_dipole_out, dipoles):
    with open(born_dipole_out, 'w') as f:
        f.write('# TimeStep dipole moment  x, y, z (|e|)\n')
        for i, dipole in enumerate(dipoles):
            f.write('{} {:.5f} {:.5f} {:.5f}\n'.format(i, *dipole))


def get_dipole_born(in_file='in.yaml'):
    with open(in_file, 'r') as stream:
        try:
            param = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)


    for key in ('dump_unfolded', 'born_file', 'born_poscar'):
        if key in param:
            param[key] = resolve_config_path(in_file, param[key])
    cfgs = read_sil(param['dump_unfolded'],
                    atom_symbols=param['atom_symbols'])
    print('cfgs were read.')

    born_dipole_out = resolve_config_path(
        in_file, param.get('born_dipole_out')
    )

    Z = get_charge_tensors(param, cfgs[0])  # shape (N_atoms, 3, 3)

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

    if born_dipole_out is not None:
        write_dipoles(born_dipole_out, dipoles)
    return dipoles


if __name__ == '__main__':
    get_dipole_born()
