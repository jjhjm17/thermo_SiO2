#!/usr/bin/env python3
"""This script prints messages for loading modules for the mlip binary."""
import sys
from parameters import server, mlip_version


def load_modules():
    """This function prints messages for loading modules."""
    if server == 'justus' and mlip_version == 3:
        print("In Justus cluster, please run the following commands in the bash "
              "by copying them: \n"
              "module purge \n"
              "module load compiler/intel/19.0  numlib/mkl/2019  mpi/impi/2019.5"
              " WARNING: the collected energy is not E0 (sigma->0), but 'F'."
              "For now, don't use mlip 3."
              )
    elif server == 'justus' and mlip_version == 2:
        print("In Justus cluster, please run the following commands in the bash "
              "by copying them: \n"
              "module purge \n"
              "module load compiler/intel  numlib/mkl  mpi/impi"
              )
    else:
        print("Error: unsupported option.")
        sys.exit()


if __name__ == '__main__':
    load_modules()
