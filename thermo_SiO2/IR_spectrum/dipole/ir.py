import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq, ifft
from scipy.signal.windows import blackman, hann
from scipy.signal import correlate
import parameter as param


def plot_debug(Y, label=''):
    plt.figure(figsize=(7,5))
    plt.plot(Y, label=label)
    plt.legend()
    plt.show()

def get_ir_dipole():
    """This function obtains the IR intensity using the autocorrelation of dipole moments, or FFT of dipole
    moments."""
    # =========================
    # User parameters
    # =========================
    # dipole_file = "dipole.dat"
    # dt_fs = 0.5                 # time step in femtoseconds
    # output_plot = "IR_spectrum.png"
    # output_data = "IR_spectrum.dat"
    dipole_file = param.dipole_file
    dt_fs = param.dt_fs
    output_plot = param.output_plot
    output_data = param.output_data
    use_window = True
    # use_window = False
    use_gradient = True  # d mu / d t
    # use_gradient = False
    if hasattr(param, 'freq_range'):
        freq_range = param.freq_range
    else:
        freq_range = None

    if hasattr(param, 'get_autocorr'):
        get_autocorr = param.get_autocorr  # get the autocorrelation function
    else:
        get_autocorr = False
    if hasattr(param, 'cut_autocorr'):
        cut_autocorr = param.cut_autocorr  # get the autocorrelation function
    else:
        cut_autocorr = 1000  # fs  # see MACE4IR paper

    # =========================
    # Load dipole data
    # =========================
    data = np.loadtxt(dipole_file)[:-1]
    # [:-1]: the length is a multile of 2, probably FFT is slightly faster?

    # time_fs = data[:, 1]
    mux = data[:, 1]
    muy = data[:, 2]
    muz = data[:, 3]

    # ref https://docs.scipy.org/doc/scipy/tutorial/fft.html
    N_sam = len(data)  # number of sample points
    ir_intensity = np.zeros(N_sam//2)  # ir_intensity: debug only
    # if get_autocorr:
    autocorr = np.zeros(N_sam, dtype=np.complex128)
    autocorr_scipy = np.zeros(N_sam * 2 - 1)

    for mu in (mux, muy, muz):
    # for mu in [muz]:

        # Remove average (important!)
        if use_gradient:
            mu = np.gradient(mu, dt_fs)
        mu -= np.mean(mu)

        # 1 cm-1 = 2.99793e10 Hz  # https://wild.life.nctu.edu.tw/class/common/energy-unit-conv-table.html
        # 1 fs = 1e-15 s = 1e-15 1/Hz = 1e-15 / (1/2.99793e10 cm-1 ) = 1e-15 * 2.99793e10 cm = 2.99793e-5
        dt_cm = dt_fs * 2.99793e-5  # unit of cm
        if use_window:
            win = blackman(N_sam)  # window
            # mu_f = fft(mu)
            mu_f = fft(mu * win)
        else:
            mu_f = fft(mu)
        freq_cm = fftfreq(N_sam, dt_cm)[:N_sam//2]  # //: floor division, 7 // 3 = 2

        ir_intensity += np.abs(mu_f[0:N_sam//2])**2

        # if get_autocorr:
        # autocorr = ifft(mu_f)  # debug
        autocorr += ifft(np.abs(mu_f)**2)
        # autocorr_shift = autocorr[0:N_sam//2]
        corr_result = correlate(mu, mu, mode='full', method='fft')
        # It changes from 0 to a large value to 0.  .../\...
        # plot_debug(corr_result)

        # autocorr_scipy += corr_result[corr_result.size//2:]  # length N_sam
        autocorr_scipy += corr_result  # length N_sam
        # ref https://stackoverflow.com/questions/643699/how-can-i-use-numpy-correlate-to-do-autocorrelation

    # ir_scipy_all = fft(autocorr_scipy)  # debug
    # plt.figure(figsize=(7,5))
    # plt.plot(ir_scipy_all)
    # plt.show()
    # sys.exit()

    # ir_scipy = fft(autocorr_scipy)[0:N_sam//2]
    half = corr_result.size // 2  # floor div
    cut = int(cut_autocorr / dt_fs)  # unit of cut: point
    autocorr_scipy = autocorr_scipy[ (half - cut) : (half + cut) ]
    plot_debug(autocorr_scipy, 'before the window')
    if use_window:
        # debug
        # autocorr_scipy = np.ones(2 * (N_sam) - 1)
        # autocorr_scipy[:(half - cut)] = 0
        autocorr_scipy *= hann(2*cut)
        # autocorr_scipy[(half + cut):] = 0
        plot_debug(hann(2*cut), 'window')
        plot_debug(autocorr_scipy, 'after the window')
    N_sam_scipy = len(autocorr_scipy)
    ir_scipy = np.abs(fft(autocorr_scipy)[0:N_sam_scipy//2])
    # np.abs is needed because fft multiplies by exp(i omega t), which is a complex number.
    freq_cm_scipy = fftfreq(N_sam_scipy, dt_cm)[:N_sam_scipy//2]  # //: floor division, 7 // 3 = 2
    # Y changes from a finite value to 0 as x increases.
    if not use_gradient:
        ir_intensity *= freq_cm **2  # omega^2
        ir_scipy *= freq_cm[ : len(ir_scipy) ] **2
    ir_intensity /= np.max(ir_intensity)
    ir_scipy /= np.max(ir_scipy)
    # =========================
    # Save spectrum
    # =========================
    # np.savetxt(output_data,
    #            np.column_stack((freq_cm, ir_intensity)),
    #            header="Frequency(cm^-1)  Intensity (AU)")
    np.savetxt(output_data,
               np.column_stack((freq_cm_scipy, ir_scipy)),
               header="Frequency(cm^-1)  Intensity (AU)")



    # =========================
    # Plot
    # =========================
    plt.figure(figsize=(7,5))
    # plt.plot(freq_cm, ir_intensity, linewidth=1.5)
    # plt.plot(freq_cm_scipy, ir_scipy, linewidth=1.5, label='IR scipy')
    plt.plot(freq_cm_scipy, ir_scipy, linewidth=1.5)
    if freq_range is not None:
        plt.xlim(freq_range)
    plt.ylim([0,1.05])
    plt.xlabel("Wavenumber (cm$^{-1}$)")
    plt.ylabel("Intensity")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_plot, dpi=300)
    plt.show()

    if get_autocorr:
        max_imag = np.max(np.abs(autocorr.imag))
        if  max_imag > 1e-6:
            print("Autocorrelation function has imaginary part (I do not know if it's okay or not.)")
            print(f'{max_imag = }')
        else:
            autocorr = autocorr.real

        autocorr /= np.max(autocorr)
        autocorr_scipy /= np.max(autocorr_scipy)

        plt.figure(figsize=(7,5))
        N_all = len(autocorr)
        x_time = np.linspace(0, (N_all-1) * dt_fs, num=N_all)
        autocorr_scipy_p = autocorr_scipy[len(autocorr_scipy)//2:] # for plot,
        # show the half
        N_all_scipy = len(autocorr_scipy_p)
        x_time_scipy = np.linspace(0, (N_all_scipy-1) * dt_fs, num=N_all_scipy)
        plt.plot(x_time, autocorr, linewidth=1.5, label='|fft|^2')
        plt.plot(x_time_scipy, autocorr_scipy_p, linewidth=1.5, label='scipy')
        plt.xlabel("Time (fs)")
        plt.ylabel("Autocorelation")
        plt.legend()
        plt.tight_layout()
        suffix = output_plot[-4:]
        plt.savefig(f'{output_plot[:-4]}_corr{suffix}', dpi=300)
        plt.show()



if __name__ == '__main__':
    get_ir_dipole()
