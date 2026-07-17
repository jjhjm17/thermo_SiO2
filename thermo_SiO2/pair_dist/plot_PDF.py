#!/usr/bin/env python3
"""This script plots the pair distribution function (PDF).  This is an edit of
An edit of Maxwell Terban's script for the calculation of pair distribution function."""
import os
import sys
from os import remove
from os.path import isfile
import numpy as np
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib import use
# from diffpy.Structure import Structure, loadStructure
# from diffpy.srreal.pdfcalculator import PDFCalculator, DebyePDFCalculator
# from diffpy.srfit.pdf.characteristicfunctions import sphericalCF
# 
# from diffpy.srreal.pdfcalculator import fftftog, fftgtof
from ase.io import read, write
from ..util.util import shell
from thermo_SiO2.mlip.read import read_SiO2
# sys.path.insert(0, '../a.total')
from a_parameters import exp_g_r, has_hydrogen
import a_parameters as param


# format = "cfg"
# config_file = "../../../../c.IR_spectrum/c.Al_SiO2/a.Al_Si_1_6/d.forceConst/a.400cfg_justus_OMP_8/sample_0/alat_Ang_14.1730044919009455/config_step_00600.cfg"
# has_hydrogen = False

def calc_PDF():
    PDF_partials = []
    PDF_totals = []
    for cfg_file in param.config_files:
        print('cfg_file: {}'.format(cfg_file))
        if cfg_file.endswith('.cfg'):
            format = 'cfg'
            kwarg = {}
        # elif cfg_file.endswith('.dump'):
        #     format = 'lammps-dump-text'
        #     specorder = ['Si', 'O', 'Al']
        #     kwarg = {'specorder': specorder}
        # elif cfg_file.endswith('.dataf'):
        #     format = 'lammps-data'
        #     specorder = param.specorder
        #     kwarg = {'specorder': specorder}
        else:
            config_ase = read_SiO2(cfg_file, atom_symbols=param.atom_symbols)

        # if hasattr(a_parameters, 'specorder'):
        #     kwarg = {'specorder': a_parameters.specorder}
        # else:
        #     kwarg = {}

        # print('Warning: for now, only the 1st configuration is used.')
        # config_file = param.config_files[0]
        # config_ase = read(config_file, index="-1", format=format, **kwarg)
        write('config.cif', config_ase, format="cif")
        shell("""
            conda run -n diffpy  python -u -m thermo_SiO2.pair_dist.PDF_py2 | tee b_plot_PDF_py2.out
            """)
        remove('config.cif')
        PDF_partials.append(np.loadtxt('PDF_partial.txt'))
        PDF_totals.append(np.loadtxt('PDF_total.txt'))
    print('Calculation of average of PDFs. ')
    PDF_partial_avg = np.average(PDF_partials, axis=0)
    PDF_total_avg = np.average(PDF_totals, axis=0)
    np.savetxt('PDF_partial.txt', PDF_partial_avg, fmt='%0.6g', header='r_l1,  gr_l1,  gr_l2,  ..., gr_l6' )
    np.savetxt('PDF_total.txt', PDF_total_avg, fmt='%0.6g', header='r_lorch_ebye, gr_lorch_debye')
    if len(cfg_file) > 1:
        PDF_total_std = np.std(PDF_totals, axis=0)
        np.savetxt('PDF_total_std.txt', PDF_total_std, fmt='%0.6g', header='r_lorch_ebye, gr_lorch_debye')


def add_subfigure_label(ax_in, string):
    ax_in.text(0.025, 0.95, string, ha='left', va='top',
               transform=ax_in.transAxes)


# def get_coord_peak(xmin, xmax, xs, ys):
#     """This function finds peak coordinate for labeling.
#     """
#     x_indices = np.where((xmin < xs) and (xs < xmax))
#     ys = ys[x_indices]
#     peaks, _ = find_peaks(ys

def add_arrow_text(ax, dx, dy, x0, y0, string, pad=-2, **kwargs):
    ax.annotate(string, xy=(x0, y0), xytext=(x0+dx, y0+dy),
                arrowprops=dict(arrowstyle="->"),
                bbox=dict(pad=pad, facecolor="none", edgecolor="none"),
                **kwargs)


def plot_PDF():

    use('pdf')
    std_avail = False
    if hasattr(param, 'PDF_folders'):
        print('Calculation of average of PDFs. ')
        std_avail = True
        if not isfile('PDF_total_std.txt'):
            print('New calculation of PDF average. ')
            PDF_partials = []
            PDF_totals = []
            for folder in param.PDF_folders:
                print('folder: {}'.format(folder))
                # PDF_partials.append(np.loadtxt(f'{folder}/PDF_partial.txt'))
                PDF_partials.append(np.loadtxt('{}/PDF_partial.txt'.format(folder)))
                # PDF_totals.append(np.loadtxt(f'{folder}/PDF_total.txt'))
                PDF_totals.append(np.loadtxt('{}/PDF_total.txt'.format(folder)))
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

    if os.path.isfile('PDF_total_std.txt'):
        std_avail = True
        PDF_total_std = np.loadtxt('PDF_total_std.txt')
        __, gr_lorch_debye_std  = PDF_total_std.T

    # config = read(param.PDF_folders[0] + '/config.cif')
    # rho_0 = len(config) / config.get_volume()

    # offset = -2
    offset = 0
    # r_baseline = np.linspace(0, 4, 2)
    # g_baseline = -4 * np.pi * rho_0 * r_baseline + offset

    size_x = 9 / 2.54
    # size_y = 7 / 2.54 * 1.0
    size_y = 7 / 2.54 * 1.7
    plt.style.use("~/program/stylelib/paper.mplstyle.no_tex")
    # fig, ax = plt.subplots(1, 1, figsize=(size_x, size_y), sharex=True)

    ylim_0 = [-1.99, 6]  # axis 0
    y_shift = -1  # shift between partial PDFs
    y0 = 0  # start of the shift of partial PDF
    num_curves = 6  # number of partial PDFs
    ylim_1 = [ylim_0[0] + (num_curves - 1) * y_shift + 1, ylim_0[1]]   # axis 1
    height_ratios = [ylim_0[1] - ylim_0[0], ylim_1[1] - ylim_1[0]]
    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(size_x, size_y), sharex=True,
                             height_ratios=height_ratios)


    ax = axes[0]
    ax.plot(r_lorch_debye, gr_lorch_debye + offset, label="Total", color='black')
    if std_avail:
        ax.fill_between(r_lorch_debye, gr_lorch_debye + gr_lorch_debye_std + offset,
                        gr_lorch_debye - gr_lorch_debye_std + offset, color='0.8')  # gray
                        # gr_lorch_debye - gr_lorch_debye_std + offset, color='0.7')  # gray
    ax.plot(dat[0], dat[1]*8 + offset, color="0.5", label="Exp")
    # ax.plot(r_baseline, g_baseline, color="0.5", linestyle='--')
    ax.legend(loc='upper right')  # matplotlib version 3.5.3
    ax.set_ylabel("$G$ (Å$^{-2}$)", labelpad=1.0)
    ax.set_ylim(ylim_0)
    ax.yaxis.set_major_locator(MultipleLocator(2))
    add_subfigure_label(ax, '(a)')


    y_shifts = np.linspace(start=y0, stop=y0+(num_curves-1)*y_shift, num=num_curves)
    ax = axes[1]
    # ax.plot(r_l3, gr_l3, label="Si-O")
    ax.plot(r_l1, gr_l3 + y_shifts[0], label=r"Si$-$O", color='C0')
    # label=r'$\mathrm{--}$' gives the same dash as $-$.
    # ax.plot(r_l5, gr_l5, label="Al-O")
    ax.plot(r_l1, gr_l5 + y_shifts[1], label="Al$-$O", linestyle='--', color='C1')
    # ax.plot(r_l6, gr_l6, label="O-O")
    # (0, (1, 5)): 'dashed', (offset, (on_off_seq))
    ax.plot(r_l1, gr_l6 + y_shifts[2], label="O$-$O", linestyle=(0, (1, 1)), color='C5')
    ax.plot(r_l1, 10*gr_l4 + y_shifts[3], label=r"Al$-$Al ($\times$10)", linestyle='--', color='C3')
    ax.plot(r_l1, gr_l1 + y_shifts[4], label="Si$-$Si", linestyle='-.', color='C4')
    # ax.plot(r_l2, gr_l2, label="Si-Al")
    ax.plot(r_l1, gr_l2 + y_shifts[5], label="Si$-$Al", linestyle='-', color='C2')
    # ax.plot(r_l4, gr_l4, label="Al-Al")
    # ax.plot(r_l4, 10*gr_l4, label=r"Al-Al ($\times$10)")
    ax.set_ylim(ylim_1)
    if has_hydrogen:
        ax.plot(r_l7, 10*gr_l7, label="H$-$O (x10)")


    ax.set_xlim(0.0,10.0)
    ax.set_xlabel("$r$ (Å)", labelpad=0.8)
    ax.set_ylabel("$G$ (Å$^{-2}$)", labelpad=0.9)
    ax.legend(loc='upper right', ncol=2)  # matplotlib version 3.5.3
    ax.yaxis.set_major_locator(MultipleLocator(2))
    add_subfigure_label(ax, '(b)')

    ax.text(1.7, 5.2, r'Si$-$O')
    # ax.annotate(r'Al$-$O', xy=(1.8, -0.6), xytext=(1.8, 1),
    #             arrowprops=dict(arrowstyle="->"),
    #             bbox=dict(pad=-2, facecolor="none", edgecolor="none"))
    ax.text(2.3, -1.2, r'Al$-$O', va='top', ha='right')
    # dx = 0.0
    # dy = 1.0
    # x0, y0 = 2.65, 1.5
    # string = r'O$-$O'
    # add_arrow_text(ax, dx, dy, x0, y0, string)
    # ax.text(3.3, -1.4, r'O$-$O', va='top', ha='center')
    ax.text(2.3, -1.9, r'O$-$O', va='top', ha='right')

    # ax.text(3.0, 1.3, r'Si$-$Si', va='bottom', ha='left')
    # add_arrow_text(ax, dx=0, dy=-1.5, x0=2.6, y0=0.45, string=r"Al$-$Al ($\times$10)",
    #                pad=-1, va='top', ha='center')
    ax.text(2.5, -3.2, r"Al$-$Al ($\times$10)", va='top', ha='right')
    ax.text(3.3, -4, r"Si$-$Si", va='bottom', ha='left')

    ax.text(4.0, 1.0, r'Si$-$O', va='bottom', ha='left')
    ax.text(3.3, -5.2, r'Si$-$Al', va='top', ha='left')
    # add_arrow_text(ax, dx=0, dy=-1.35, x0=4.9, y0=0.3, string=r"Si$-$Si",
    #                pad=-1, va='top', ha='center')

    fig.subplots_adjust(
        top=0.985, bottom=0.060, left=0.085, right=0.985, hspace=0.0,
        wspace=0.2
    )
    # fig.savefig('fig__pair_dist_no_exp.pdf', dpi=1000, transparent=True)
    fig.align_ylabels()
    fig.savefig('fig__pair_dist.pdf', dpi=1000, transparent=True)

    # fig.savefig('fig_pair_dist.pdf', transparent=True)
    # if make_detail_figures:
    # plt.show()



if __name__ == '__main__':
    plot_PDF()
