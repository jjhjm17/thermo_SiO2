#!/usr/bin/env python3
"""This function obtaions the radial distrubiton function of O. The O-O
peak has effect on the diffusion coefficient."""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from ase.io import read
from ase.geometry.analysis import Analysis
# from ase.atoms import Atoms
from ase.atom import Atom
from ..mlip.read import read_SiO2
from a_parameters import file, rmax, nbins, index, atom_symbols, elements_rdf


def get_rdf():
    """This function obtains the radial distribution function."""
    # configs = read(file, index=index, format='lammps-dump-text',
    #             specorder=specorder)
    # if isinstance(configs[0], Atom):
    #     configs= [configs]  # configs
    configs = read_SiO2(file, atom_symbols=atom_symbols, index=index)
    print(f'{len(configs)} configs., {configs[0] = }')
    ana = Analysis(configs)

    rdf_X_data = ana.get_rdf(rmax=rmax, nbins=nbins, return_dists=True,
                             elements=elements_rdf)
    # print(f'{rdf_X_data = }')
    rdf_X_data = np.array(rdf_X_data)
    # print(f'{rdf_X_data[:, 0] = }')
    rdf = np.average(rdf_X_data[:, 0], axis=0)
    X = rdf_X_data[0][1]
    # print(f'{rdf = }')

    # rdf = rdf_X[0]
    # X = rdf_X[1]
    # print(f'{rdf = }')

    fig, ax = plt.subplots()
    ax.plot(X, rdf, color='green')
    ax.set_xlim([2, rmax])
    ax.set_ylim([0, 4])
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.set_xlabel(r'$r$ (Ang)')
    ax.set_ylabel(r'$g_\mathrm{OO}$ (Ang)')
    fig.savefig('fig-rdf.pdf')
    plt.show()


if __name__ == '__main__':
    get_rdf()
