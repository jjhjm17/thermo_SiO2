"""This script contains functions for writing SiO2 configurations.

Be carefule:
    the mlippy.atms.ase_savecfgs is modified from the mlippy in mlip-2.
    mlippy/cfgs.py file is modified by J. H. Jung.
    so  that mlippy  ase_loadcfgs can read cfgs files written by the function  ase_savecfgs.

"""
from thermo_SiO2._vendor.mlippy.atms import ase_savecfgs
from ..util.SiO2_parameter import atomic_number_to_POTCAR_order, Si_O_Al_atomic_number_to_POTCAR_order
from thermo_SiO2.util.SiO2_parameter import Atom_order


def write_cfg_SiO2(file, configs, atom_symbols='Si O H Al', desc=None):
    """This function writes a mlip cfg file of SiO2 structure and changes the
    atom type correctly. A list of Atoms-type objects is returned."""
    # 'configs' is a list of Atoms-type objects.

    if type(configs).__name__ == 'Atoms':
        configs = [configs]
        # If there is one configuration, convert to a list of configs.
    # if atom_symbols == 'Si O H Al': 
    #     to_POTCAR_order = atomic_number_to_POTCAR_order
    # elif atom_symbols == 'Si O Al':
    #     to_POTCAR_order = Si_O_Al_atomic_number_to_POTCAR_order
    # else:
    #     print('Error: please code more for this atom_symbols.')
    atom_order = Atom_order(atom_symbols)
    to_POTCAR_order = atom_order.atomic_number_to_POTCAR_order()

    for config in configs:
        atomic_numbers = config.get_atomic_numbers()
        atomic_numbers = [to_POTCAR_order[number] for number in
                          atomic_numbers]
        config.set_atomic_numbers(atomic_numbers)

    ase_savecfgs(file, configs, desc=desc)
    return configs
