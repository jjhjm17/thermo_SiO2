#!/usr/bin/env python3
"""This script plots the pair distribution function (PDF).  This is an edit of
An edit of Maxwell Terban's script for the calculation of pair distribution function.
Because the default channel of anaconda is now commercial, we use the conda-forge channel, and so we can use only python 2.7. We cannot use python 3."""
from __future__ import print_function, division

import numpy as np
import matplotlib.pyplot as plt
from diffpy.Structure import Structure, loadStructure
from diffpy.srreal.pdfcalculator import PDFCalculator, DebyePDFCalculator
from diffpy.srfit.pdf.characteristicfunctions import sphericalCF

from diffpy.srreal.pdfcalculator import fftftog, fftgtof
from a_parameters import has_hydrogen


def Lorch(Qmax, Q):
    """
    Lorch function. It minimizes the truncation errors associated with Qmax.
    """
    delta = np.pi/Qmax
    return np.sin(Q*delta)/(Q*delta)


def Lorch_to_g(ftrunc, qmax, pad):
    """
    Fourier transform F(Q) data to G(r) data with the Lorch function applied
    """
    
    qstep = 0.01
    
    newq = np.interp(np.arange(0.0,qmax,0.01), ftrunc[0], ftrunc[1])
    q = np.arange(0.0,qmax,qstep)
    fq = newq*Lorch(qmax, q)
    
    q = np.pad(q, (0,pad), 'constant')
    fq = np.pad(fq, (0,pad), 'constant')
    
    gr, rstep = fftftog(fq, qstep)
    r = np.arange(len(gr)) * rstep
    
    return r, gr


def g_to_f(gtrunc, rmax, pad):
    """
    Fourier transform F(Q) data to G(r) data with the Lorch function applied
    """
    
    rstep = 0.01
    
    newr = np.interp(np.arange(0.0,rmax,0.01), gtrunc[0], gtrunc[1])
    r = np.arange(0.0,rmax,rstep)
    gr = newr
    
    r = np.pad(r, (0,pad), 'constant')
    gr = np.pad(gr, (0,pad), 'constant')
    
    fq, qstep = fftgtof(gr, rstep)
    q = np.arange(len(fq)) * qstep
    
    return q, fq


def calc_PDF_py2():
    """This function reads config.cif and writes PDF_partial.txt, PDF_total.txt. python 2 is used."""

    struct_file = 'config.cif'

    # struc = loadStructure('Al_SiO2.cif')
    struc = loadStructure(struct_file)
    struc.Uisoequiv = 0.001

    formula = ["Si","Al","O"]
    #biso = [0.001,0.001,0.001]

    #######Calculate using real-space calculator, i.e. form factors equal f(0)#########
    cfg = { 'qmax' : 0,'qmin': 0,'rmin' : 0.01,'rmax' : 20.0,'qdamp': 0.01,'qbroad': 0.01}

    pc1 = PDFCalculator(**cfg) 
    r1, g1 = pc1(struc)

    #######Calculate using Debye calculator, i.e. form factors equal f(Q)#########
    cfg2 = { 'qmax' : 50,'qmin': 0,'rmin' : 0.01,'rmax' : 20.0,'qdamp': 0.01,'qbroad': 0.01}

    pc2 = DebyePDFCalculator(**cfg2) 
    r2, g2 = pc2(struc)

    ###########Approximate small-angle scattering contribution#############
    struc.Uisoequiv = 1.0

    pc2 = DebyePDFCalculator(**cfg2) 
    r3, g3 = pc2(struc)

    plt.figure(0)
    plt.plot(r1,(g1)*sphericalCF(r1,15.0), label="PDFcalculator")
    plt.plot(r2,(g2-g3)*sphericalCF(r2,15.0), label="DebyePDFcalculator")
    plt.plot(r2,g2, label="DebyePDFcalculator g2")
    plt.plot(r2,sphericalCF(r2,15.0), label="sphericalCF")
    # print(f'r3 = {r3}')
    # print(f'g3 = {g3}')
    plt.plot(r3,g3, label='small angle scattering')

    plt.xlim(0.0,10.0)
    plt.legend()

    plt.savefig('fig__a_small_angle_scatttering.pdf', transparent=True)
    # plt.show()

    ########Fourier transform the calculated PDFs#############
    plt.figure(1)
    q1, fq1 = g_to_f([r1,(g1)*sphericalCF(r1,15.0)], 100.0, 1000)
    q2, fq2 = g_to_f([r2,(g2-g3)*sphericalCF(r2,15.0)], 100.0, 1000)
    plt.plot(q1,fq1, label="PDFcalculator")
    plt.plot(q2,fq2, label="DebyePDFcalculator")

    plt.xlim(0.0,50.0)
    plt.ylim(-5,7)
    plt.legend()
    plt.savefig('fig__b_DebyePDFcalc.pdf', transparent=True)
    # plt.show()

    ########Fourier transform the calculated reduced structure functions with Lorch function applied#############

    plt.figure(2)
    r_l1, gr_l1 = Lorch_to_g([q1,fq1], 30.0, 10000)
    r_l2, gr_l2 = Lorch_to_g([q2,fq2], 30.0, 10000)
    r_lorch_debye, gr_lorch_debye = r_l2, gr_l2
    # dat = np.loadtxt(exp_g_r).T

    plt.plot(r_l1, gr_l1, label="PDFcalculator")
    plt.plot(r_l2, gr_l2, label="DebyePDFcalculator")
    # plt.plot(dat[0], dat[1]*8-2.0, color="black", label="experiment")

    plt.xlim(0.0,10.0)
    plt.legend()

    plt.savefig('fig__c_Fourier_transform_Lorch.pdf', transparent=True)
    # plt.show()

    plt.figure(3)
    # struc = loadStructure('Al_SiO2.cif')
    struc = loadStructure(struct_file)
    struc.Uisoequiv = 0.001

    #######Calculate using Debye calculator, i.e. form factors equal f(Q)#########
    # cfg2 = { 'qmax' : 50,'qmin': 0,'rmin' : 0.01,'rmax' : 20.0,'qdamp': 0.01,'qbroad': 0.01}
    cfg2 = { 'qmax' : 30,'qmin': 0,'rmin' : 0.01,'rmax' : 20.0,'qdamp': 0.01,'qbroad': 0.01}

    pc2 = DebyePDFCalculator(**cfg2) 
    pc2.setTypeMask('Si', 'Si', True, others=False)
    r1, g1 = pc2(struc)
    pc2.setTypeMask('Si', 'Al', True, others=False)
    r2, g2 = pc2(struc)
    pc2.setTypeMask('Si', 'O', True, others=False)
    r3, g3 = pc2(struc)
    pc2.setTypeMask('Al', 'Al', True, others=False)
    r4, g4 = pc2(struc)
    pc2.setTypeMask('Al', 'O', True, others=False)
    r5, g5 = pc2(struc)
    pc2.setTypeMask('O', 'O', True, others=False)
    r6, g6 = pc2(struc)
    if has_hydrogen:
        pc2.setTypeMask('H', 'O', True, others=False)
        r7, g7 = pc2(struc)

    ###########Approximate small-angle scattering contribution#############
    struc.Uisoequiv = 10.0

    pc2 = DebyePDFCalculator(**cfg2) 
    pc2.setTypeMask('Si', 'Si', True, others=False)
    r12, g12 = pc2(struc)
    pc2.setTypeMask('Si', 'Al', True, others=False)
    r22, g22 = pc2(struc)
    pc2.setTypeMask('Si', 'O', True, others=False)
    r32, g32 = pc2(struc)
    pc2.setTypeMask('Al', 'Al', True, others=False)
    r42, g42 = pc2(struc)
    pc2.setTypeMask('Al', 'O', True, others=False)
    r52, g52 = pc2(struc)
    pc2.setTypeMask('O', 'O', True, others=False)
    r62, g62 = pc2(struc)
    if has_hydrogen:
        pc2.setTypeMask('H', 'O', True, others=False)
        r72, g72 = pc2(struc)

    plt.plot(r1,(g1-g12)*sphericalCF(r2,15.0), label="Si-Si")
    plt.plot(r2,(g2-g22)*sphericalCF(r2,15.0), label="Si-Al")
    plt.plot(r3,(g3-g32)*sphericalCF(r2,15.0), label="Si-O")
    plt.plot(r4,(g4-g42)*sphericalCF(r2,15.0), label="Al-Al")
    plt.plot(r5,(g5-g52)*sphericalCF(r2,15.0), label="Al-O")
    plt.plot(r6,(g6-g62)*sphericalCF(r2,15.0), label="O-O")
    if has_hydrogen:
        plt.plot(r7,(g7-g72)*sphericalCF(r2,15.0), label="H-O")

    plt.xlim(0.0,10.0)
    plt.legend()

    plt.savefig('fig__d_approximate_small_angle.pdf', transparent=True)
    # plt.show()

    ########Fourier transform the calculated PDFs to F(Q)#############

    q1, fq1 = g_to_f([r1,(g1-g12)*sphericalCF(r2,15.0)], 100.0, 1000)
    q2, fq2 = g_to_f([r2,(g2-g22)*sphericalCF(r2,15.0)], 100.0, 1000)
    q3, fq3 = g_to_f([r3,(g3-g32)*sphericalCF(r2,15.0)], 100.0, 1000)
    q4, fq4 = g_to_f([r4,(g4-g42)*sphericalCF(r2,15.0)], 100.0, 1000)
    q5, fq5 = g_to_f([r5,(g5-g52)*sphericalCF(r2,15.0)], 100.0, 1000)
    q6, fq6 = g_to_f([r6,(g6-g62)*sphericalCF(r2,15.0)], 100.0, 1000)
    if has_hydrogen:
        q7, fq7 = g_to_f([r7,(g7-g72)*sphericalCF(r2,15.0)], 100.0, 1000)

    ########Fourier transform the calculated reduced structure functions with Lorch function applied#############

    r_l1, gr_l1 = Lorch_to_g([q1, fq1], 30.0, 10000)
    r_l2, gr_l2 = Lorch_to_g([q2, fq2], 30.0, 10000)
    r_l3, gr_l3 = Lorch_to_g([q3, fq3], 30.0, 10000)
    r_l4, gr_l4 = Lorch_to_g([q4, fq4], 30.0, 10000)
    r_l5, gr_l5 = Lorch_to_g([q5, fq5], 30.0, 10000)
    r_l6, gr_l6 = Lorch_to_g([q6, fq6], 30.0, 10000)
    if has_hydrogen:
        r_l7, gr_l7 = Lorch_to_g([q7, fq7], 30.0, 10000)
    PDF_partial = np.column_stack((r_l1,  gr_l1, gr_l2, gr_l3, gr_l4, gr_l5,
                                   gr_l6))
    PDF_total = np.column_stack((r_lorch_debye, gr_lorch_debye)) 

    for r_l in [r_l2, r_l3, r_l4, r_l5, r_l6]:
        if not np.allclose(r_l, r_l1):
            print('Error: r_l is different.')
            sys.exit()

    np.savetxt('PDF_partial.txt', np.column_stack((r_l1,  gr_l1, gr_l2, gr_l3,  gr_l4, gr_l5, gr_l6)), fmt='%0.6g', header='r_l1,  gr_l1,  gr_l2,  ..., gr_l6' )
    np.savetxt('PDF_total.txt', np.column_stack((r_lorch_debye, gr_lorch_debye)), fmt='%0.6g', header='r_lorch_ebye, gr_lorch_debye')

    return PDF_partial, PDF_total


if __name__ == '__main__':
    calc_PDF_py2()
