#!/bin/bash
# python3 -u -m thermo_SiO2.active_learn.lammps.create_folders | tee ./b.create_folders.out
python3 -u  create_folders.py | tee ./b.create_folders.out
exit
