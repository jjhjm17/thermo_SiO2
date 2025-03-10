#!/usr/bin/env python3
"""This script prints messages for loading modules for the mlip binary."""
import sys
from parameters import server, mlip_version


def load_modules():
    """This function prints messages for loading modules."""
    print(f'{server = }, {mlip_version = }')
    if server == 'justus' and mlip_version == 3:
        print("In Justus cluster, please run the following commands in the bash "
              "by copying it: \n"
              "load_module_mlp3_lammps \n"
              " WARNING: the collected energy is not E0 (sigma->0), but 'F'."
              "For now, don't use mlip 3."
              )
    elif server == 'justus' and mlip_version == 2:
        print("In Justus cluster, please run the following commands in the bash "
              "by copying them: \n"
              "module purge \n"
              "module load compiler/intel  numlib/mkl  mpi/impi"
              )
    elif server == 'fritz' and mlip_version == 2:
        print("In Fritz cluster, please run the following commands in the bash "
              "by copying them: \n"
              "module purge \n"
              "module load intel mkl intelmpi \n"
              "module load tbb  # intel Threading Building Blocks"
              )
    elif server == 'fritz' and mlip_version == 3:
        print("In Fritz cluster, please run the following commands in the bash "
              "by copying them: \n"
              "load_module_mlp3_lammps \n"
              " WARNING: the collected energy is not E0 (sigma->0), but 'F'.")
    else:
        print("Error: unsupported option.")
        sys.exit()


if __name__ == '__main__':
    load_modules()
