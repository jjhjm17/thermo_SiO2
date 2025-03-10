#!/usr/bin/env python3
"""This script plots the pair distribution function (PDF).  This is an edit of
An edit of Maxwell Terban's script for the calculation of pair distribution function."""
from os.path import isfile
# import sys
import numpy as np
import matplotlib.pyplot as plt

# from diffpy.Structure import Structure, loadStructure
# from diffpy.srreal.pdfcalculator import PDFCalculator, DebyePDFCalculator
# from diffpy.srfit.pdf.characteristicfunctions import sphericalCF
# 
# from diffpy.srreal.pdfcalculator import fftftog, fftgtof
from ase.io import read, write
from ..util.util import shell
# sys.path.insert(0, '../a.total')
from a_parameters import exp_g_r, has_hydrogen
import a_parameters as param


# format = "cfg"
# config_file = "../../../../c.IR_spectrum/c.Al_SiO2/a.Al_Si_1_6/d.forceConst/a.400cfg_justus_OMP_8/sample_0/alat_Ang_14.1730044919009455/config_step_00600.cfg"
# has_hydrogen = False

def calc_PDF():
    if param.config_files[0].endswith('.cfg'):
        format = 'cfg'
        kwarg = {}
    elif param.config_files[0].endswith('.dump'):
        format = 'lammps-dump-text'
        specorder = ['Si', 'O', 'Al']
        kwarg = {'specorder': specorder}

    # if hasattr(a_parameters, 'specorder'):
    #     kwarg = {'specorder': a_parameters.specorder}
    # else:
    #     kwarg = {}

    print('Warning: for now, only the 1st configuration is used.')
    config_file = param.config_files[0]
    config_ase = read(config_file, index="-1", format=format, **kwarg)
    write('config.cif', config_ase, format="cif")

    shell("""
          conda run -n diffpy  python -u -m thermo_SiO2.pair_dist.PDF_py2 | tee b_plot_PDF_py2.out
          """)


def main():

    std_avail = False
    if hasattr(param, 'PDF_folders'):
        print('Calculation of average of PDFs. ')
        std_avail = True
        if not isfile('PDF_total_std.txt'):
            print('New calculation of PDF average. ')
            PDF_partials = []
            PDF_totals = []
            for folder in param.PDF_folders:
                print(f'folder: {folder}')
                PDF_partials.append(np.loadtxt(f'{folder}/PDF_partial.txt'))
                PDF_totals.append(np.loadtxt(f'{folder}/PDF_total.txt'))
            PDF_partial_avg = np.average(PDF_partials, axis=0)
            PDF_total_avg = np.average(PDF_totals, axis=0)
            PDF_total_std = np.std(PDF_totals, axis=0)
            np.savetxt('PDF_partial.txt', PDF_partial_avg, fmt='%0.6g', header='r_l1,  gr_l1,  gr_l2,  ..., gr_l6' )
            np.savetxt('PDF_total.txt', PDF_total_avg, fmt='%0.6g', header='r_lorch_ebye, gr_lorch_debye')
            np.savetxt('PDF_total_std.txt', PDF_total_std, fmt='%0.6g', header='r_lorch_ebye, gr_lorch_debye')

    if not isfile('PDF_total.txt'):
        print('New calculation of PDF. ')

        # PDF_partial, PDF_total = calc_PDF()
        calc_PDF()
    else:
        print('Previous calculation exists and we read PDF_total.txt and PDF_partial.txt. ')

    PDF_partial = np.loadtxt('PDF_partial.txt')
    PDF_total = np.loadtxt('PDF_total.txt')

    dat = np.loadtxt(exp_g_r).T

    r_l1,  gr_l1, gr_l2, gr_l3, gr_l4, gr_l5, gr_l6 = PDF_partial.T
    r_lorch_debye, gr_lorch_debye  = PDF_total.T
    if std_avail:
        PDF_total_std = np.loadtxt('PDF_total_std.txt')
        __, gr_lorch_debye_std  = PDF_total_std.T

    config = read(param.PDF_folders[0] + '/config.cif')
    # rho_0 = len(config) / config.get_volume()

    offset = -2
    # r_baseline = np.linspace(0, 4, 2)
    # g_baseline = -4 * np.pi * rho_0 * r_baseline + offset

    size_x = 9 / 2.54
    # size_y = 7 / 2.54 * 1.0
    size_y = 7 / 2.54 * 1.4
    plt.style.use("~/program/stylelib/paper.mplstyle")
    # fig, ax = plt.subplots(1, 1, figsize=(size_x, size_y), sharex=True)
    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(size_x, size_y), sharex=True)


    ax = axes[0]
    ax.plot(r_lorch_debye, gr_lorch_debye + offset, label="Total", color='black')
    if std_avail:
        ax.fill_between(r_lorch_debye, gr_lorch_debye + gr_lorch_debye_std + offset,
                        gr_lorch_debye - gr_lorch_debye_std + offset, color='0.8')  # gray
                        # gr_lorch_debye - gr_lorch_debye_std + offset, color='0.7')  # gray
    ax.plot(dat[0], dat[1]*8 + offset, color="0.5", label="Exp")
    # ax.plot(r_baseline, g_baseline, color="0.5", linestyle='--')


    ax = axes[1]
    # ax.plot(r_l3, gr_l3, label="Si-O")
    ax.plot(r_l1, gr_l3, label=r"Si$-$O")
    # label=r'$\mathrm{--}$' gives the same dash as $-$.
    # ax.plot(r_l5, gr_l5, label="Al-O")
    ax.plot(r_l1, gr_l5, label="Al$-$O")
    # ax.plot(r_l6, gr_l6, label="O-O")
    ax.plot(r_l1, gr_l6, label="O$-$O")
    ax.plot(r_l1, gr_l1, label="Si$-$Si")
    # ax.plot(r_l2, gr_l2, label="Si-Al")
    ax.plot(r_l1, gr_l2, label="Si$-$Al")
    # ax.plot(r_l4, gr_l4, label="Al-Al")
    # ax.plot(r_l4, 10*gr_l4, label=r"Al-Al ($\times$10)")
    ax.plot(r_l1, 10*gr_l4, label=r"Al$-$Al ($\times$10)")
    if has_hydrogen:
        ax.plot(r_l7, 10*gr_l7, label="H$-$O (x10)")


    ax.set_xlim(0.0,10.0)
    ax.set_ylim([-3.5, 6])
    ax.set_xlabel("$r$ (Å)", labelpad=1.0)
    ax.set_ylabel("$G$ (Å$^{-2}$)", labelpad=1.0)
    ax.legend(loc='upper right', ncol=2)  # matplotlib version 3.5.3
    fig.subplots_adjust(
        top=0.98, bottom=0.10, left=0.09, right=0.985, hspace=0.2, wspace=0.2
    )

    # fig.savefig('fig__pair_dist_no_exp.png', dpi=1000, transparent=True)
    # fig.savefig('fig__pair_dist_no_exp.pdf', dpi=1000, transparent=True)
    fig.savefig('fig__pair_dist.pdf', dpi=1000, transparent=True)

    # fig.savefig('fig_pair_dist.pdf', transparent=True)
    # if make_detail_figures:
    plt.show()


if __name__ == '__main__':
    main()
