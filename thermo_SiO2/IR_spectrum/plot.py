import numpy as np


def loadtxt_exp(fname, factor=1):
    """This function loads experimental IR."""
    data = np.loadtxt(fname)
    T = data[:,1] / 100  # Transmittance, T = I/I0, % to unitless
    # A: absorbance
    # A = log(I0/I) = A = epsilon l c
    # Beer-Lambert law
    # https://en.wikipedia.org/wiki/Beer%E2%80%93Lambert_law
    # A = log(1/T) = -log(T)
    A = -np.log10(T)
    # factor = 0.5  # arbitrary unit, for comparison with calc
    A -= np.min(A)  # baseline correction

    # normalize
    # print(f'{np.max(A[1:]) =}')  # The 1st value is 0 in a file.
    # A /= np.max(A[1:])  # normalize to 1 because density is different, probably.

    return np.column_stack((data[:,0], factor * A))


def loadtxt_ir(fname):
    """This function load IR data."""
    return np.loadtxt(fname, skiprows=1)


def get_calc_average(folder_sample_in, normalize=False, file='IR-Spectrum.dat'):
    """This function calculates average of calculation. if
    scale_factor==None, the scale factor is calculated.
    normalize==True: normalize so that the maximum value is 1"""
    spectra_all = []
    # print(f'{len(folder_sample) =}')
    # for i_sample in [1, 2, 3, 4]:
    for i_sample in range(1, len(folder_sample_in) + 1):  # 1 to len(folder_sample)
        data = loadtxt_ir(f'{folder_sample_in[i_sample]}/{file}')
        spectra_all.append(data[:,1])
    spectra_all = np.array(spectra_all)
    calc_avg = np.average(spectra_all, axis=0)
    calc_std = np.std(spectra_all, axis=0)
    X = data[:,0]
    if normalize:
        max_avg = np.max(calc_avg)
        calc_avg /= max_avg
        calc_std /= max_avg
        spectra_all /= max_avg
    return spectra_all, calc_avg, X

