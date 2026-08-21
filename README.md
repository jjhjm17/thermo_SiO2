# thermo_SiO2

Functions for the calculation of amorphous and mesoporous SiO2

## Features

- generate lammps input files for amorphous silica (Al2O3)x/2(SiO2)y
- active learning
- diffusion analysis
- pair distribution functions
- infrared spectrum analysis
- MTP utilities
- ...

## Installation

pip install ...

or

git clone ...

## Requirements

Python >= ...

ASE
NumPy
SciPy
Pandas
...

Tested versions

   Python 3.11, numpy 2.0.2, pandas 2.2.3, ase 3.23.0
   OS: AlmaLinux 8
   A high-performance computing cluster is required to run VASP, LAMMPS, and MLIP
   calculations.
   For IR: phonopy 4.4.0, phonopy-ir 28 May 2023,
       hiphive 1.4 (hiphive is needed if the effective hessians are used)
   Optional: lammps-mtp-kokkos 1.0.0 (an optimized MTP implementation for LAMMPS,
       https://github.com/RichardZJM/lammps-mtp-kokkos)


## Usage

- generate lammps input files for amorphous silica (Al2O3)x/2(SiO2)y
   see examples/melt_quench/
- visualize mlip version 2, 3 cfg files, and lammps dump files using
  ase gui with default 'Si O H Al' atomic symbols:
    visi a.cfg
    visi a.dump
    visi a.dataf
  For more information, see visi -h

## Repository structure

IR_spectrum/
pair_dist/
active_learn/
...

## Third-party software

This repository vendors a modified copy of mlippy from MLIP-2.

The modifications are limited to output formatting.

This repository contains a modified copy of mlippy from the MLIP-2 project under the BSD 2-Clause License. See _vendor/mlippy/LICENSE.mlippy .

## Citation

None 

## License

This project is licensed under the MIT License.

Third-party code under _vendor is licensed separately.
