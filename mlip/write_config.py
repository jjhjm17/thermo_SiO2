"""This script contains functions for writing SiO2 configurations."""
from mlippy.atms import ase_savecfgs
from ..util.SiO2_parameter import atomic_number_to_POTCAR_order


def write_cfg_SiO2(file, configs, desc):
    """This function writes a mlip cfg file of SiO2 structure and changes the
    atom type correctly. A list of Atoms-type objects is returned."""
    # 'configs' is a list of Atoms-type objects.
    if type(configs).__name__ == 'Atoms':
        configs = [configs]
        # If there is one configuration, convert to a list of configs.
    for config in configs:
        atomic_numbers = config.get_atomic_numbers()
        atomic_numbers = [atomic_number_to_POTCAR_order[number] for number in
                          atomic_numbers]
        config.set_atomic_numbers(atomic_numbers)

    ase_savecfgs(file, configs, desc=desc)
    return configs
