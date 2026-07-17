#!/usr/bin/env python3
"""This script calculates the radial density profile of H2O inside pore."""
import sys
import numpy as np
from numpy.linalg import norm
from scipy.ndimage import gaussian_filter1d
from sklearn.neighbors import KernelDensity
from ase.io import read
# from matplotlib import use
# import matplotlib.pyplot as plt
# from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from a_parameters import config_without_H2O_file, config_with_H2O_file, hole_radius, num_points, sigma
import a_parameters
from ..mlip.read_config import read_SiO2_dump


config_without_H2O = read(config_without_H2O_file, format='lammps-data')

if hasattr(a_parameters, 'num_configs'):
    index = f':{a_parameters.num_configs}'
else:
    index = ':'
configs = read_SiO2_dump(config_with_H2O_file, index=index)  # debug
print(f'{config_without_H2O =}')
print(f'{configs[:2] =}')
# sys.exit()

config = configs[0]
O_list = []
H_list = []
# r_s: radial distances
for id, symbol in enumerate(config.get_chemical_symbols()):
    # print(f'{id = } {symbol =} {len(struct_without_H2O) = }')
    if (id >= len(config_without_H2O) and symbol == 'O'):
        O_list.append(id)
    elif (id >= len(config_without_H2O) and symbol == 'H'):
        H_list.append(id)
print(f'{len(O_list) = }')
print(f'{len(H_list) = }')
print(f'Please check the number of H2O: {3 * len(O_list)} in log.lammps')

symbols = config.get_chemical_symbols()
with open(config_with_H2O_file, 'r') as fin:
    for ind, line in enumerate(fin):
        if ind == 5:
            xlo_bound, xhi_bound, xy = line.rstrip().split()
        elif ind == 6:
            ylo_bound, yhi_bound, xz = line.rstrip().split()
        elif ind >= 7:
            break
xlo_bound = float(xlo_bound)
xhi_bound = float(xhi_bound)
center_x = xlo_bound + (xhi_bound - xlo_bound) / 2
ylo_bound = float(ylo_bound)
yhi_bound = float(yhi_bound)
center_y = ylo_bound + (yhi_bound - ylo_bound) / 2
print(f'{center_x =}, {center_y =}')
center_xy = np.array((center_x, center_y))

O_dist = []
H_dist = []
for config in configs:
    for id, position in enumerate(config.get_positions()):
        if id in O_list:
            O_dist.append( norm( position[:2] - center_xy))
        elif id in H_list:
            H_dist.append( norm( position[:2] - center_xy))

print(f'{len(O_dist) = }')
print(f'{len(H_dist) = }')
# debug
# O_dist = [2]

x_max = 1.5 * hole_radius
dx = x_max / (num_points)
x_grid = np.linspace(dx, x_max, num=num_points)[:, np.newaxis]
# exclude 0 because later we divide by 0

height = config.get_cell()[2,2]
print(f'{height = }')



# half_width_gaussian = round(3*sigma/dx) * dx
# num_gaussian = int((2 * half_width_gaussian) / dx + 1)
# x_gaussian = np.linspace(- half_width_gaussian, half_width_gaussian, num=
#                          num_gaussian)
# print(f'{dx = }, {half_width_gaussian = }, {x_gaussian = }') 
# 
# # Ref
# # https://stackoverflow.com/questions/24148902/python-convolution-with-a-gaussian
# gaussian = 1 / (sigma * np.sqrt(2 * np.pi)) * np.exp(-0.5 * (x_gaussian / sigma) ** 2) 
# print(f'{gaussian = }')
# 
# # Normalize so that we do not change the density.
# area = np.sum(gaussian) * dx
# gaussian /= area
# print(f'{np.sum(gaussian) * dx = }')
# 
# O_dist_conv = np.convolve(O_dist

def get_radial_density(dist, x_grid, num_atoms_in_cell):
    dist = np.array(dist)
    # Ref.: https://scikit-learn.org/stable/auto_examples/neighbors/plot_kde_1d.html#sphx-glr-auto-examples-neighbors-plot-kde-1d-py
    dist = dist[:, np.newaxis]
    kde = KernelDensity(kernel='gaussian', bandwidth=sigma).fit(dist)
    log_dens = kde.score_samples(x_grid)
    dens = np.exp(log_dens)  # distribution, area=1
    # # average, radial density
    dens = dens * num_atoms_in_cell
    dens *= 1e3  # Ang-3 to nm-3
    for index, x_value in enumerate(x_grid[:, 0]):
        # if x_value < 1e-6:
        #     dens[index] = np.NaN
        # else:
        dens[index] = dens[index] / (2 * np.pi * x_value *  height + 1e-3)
    return dens


# O_dist = np.array(O_dist)
# # Ref.: https://scikit-learn.org/stable/auto_examples/neighbors/plot_kde_1d.html#sphx-glr-auto-examples-neighbors-plot-kde-1d-py
# O_dist = O_dist[:, np.newaxis]
# kde_O = KernelDensity(kernel='gaussian', bandwidth=sigma).fit(O_dist)
# log_dens_O = kde_O.score_samples(x_grid)
# dens_O = np.exp(log_dens_O)

dens_O = get_radial_density(O_dist, x_grid, num_atoms_in_cell=len(O_list))
dens_H = get_radial_density(H_dist, x_grid, num_atoms_in_cell=2*len(O_list))


np.savetxt('density.txt', np.column_stack((x_grid[:, 0], dens_O, dens_H, dens_H/2)),
           header='# X_grid (Ang), number density of O, rho of H, (rho of H) / 2 ')

# # Single column 9 cm / 2.54 (cm/inch)
# size_x = 9 / 2.54
# size_y = 6 / 2.54 * 1.0
# 
# fig, ax = plt.subplots(1,1, figsize=(size_x, size_y), sharex=True)
