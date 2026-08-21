import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq
from scipy.signal.windows import hann
from scipy.signal import correlate


def plot_debug(Y, verbose, label=''):
    if verbose:
        plt.figure(figsize=(7, 5))
        plt.plot(Y, label=label)
        plt.legend()
        plt.show()


def compute_ir_from_dipole_components(mu_components, dt_fs, use_window=True,
                                       use_gradient=True, cut_autocorr=1000,
                                       normalize=True, verbose=False):
    """Compute an isotropic IR spectrum from a set of dipole-moment (or
    dipole-like, e.g. an O-H bond vector times a formal charge) Cartesian
    components, via the dipole-autocorrelation-function (DACF) -> window
    -> FFT method.

    This is the core computation shared by get_ir_dipole() (whole-cell
    dipole moment) and get_OH_dipoles_ir.get_OH_dipoles_ir() (per-O-H-group
    partial dipole).

    Parameters
        mu_components : sequence of 1D arrays, e.g. [mux, muy, muz],
            each of length N_sam, evenly sampled every dt_fs.
        dt_fs : float
            Time step in femtoseconds between successive samples.
        use_window : bool
            Apply a Hann window to the truncated autocorrelation function
            before the final FFT, to suppress spectral leakage/ringing
            from the truncation at cut_autocorr.
        use_gradient : bool
            If True, use d(mu)/dt (the "velocity dipole") instead of mu
            itself. This avoids double-counting the omega^2 factor that
            otherwise has to be applied by hand (see below), since
            FT[dmu/dt] already carries one factor of i*omega relative to
            FT[mu].
        cut_autocorr : float
            Autocorrelation cutoff length, in fs. The frequency
            resolution of the resulting spectrum is set by this cutoff
            (~ 1 / (cut_autocorr * c)), not by the full trajectory
            length; see the MACE4IR paper. Silently capped at the
            trajectory length if a longer value is given.
        normalize : bool
            If True (default), rescale so that max(ir) == 1 -- this is
            what get_ir_dipole() has always done, and is fine when only
            a single spectrum is plotted. Set to False to keep the raw
            (unscaled) intensity, which is needed to compare relative
            intensities *between* separately-computed spectra (e.g. the
            partial spectrum of one O-H bond vs. another, in
            get_OH_dipoles_ir.py) -- with normalize=True each spectrum's
            own peak is forced to 1, so such comparisons would be
            meaningless.
        verbose : bool
            If True, show debug plots of the autocorrelation function
            before/after windowing (useful when tuning cut_autocorr).

    Returns
        freq_cm : 1D array
            Wavenumber (cm^-1), positive-frequency half only.
        ir : 1D array
            Intensity. Normalized so that max(ir) == 1 if normalize is
            True; otherwise the raw (arbitrary but internally
            consistent) intensity from the FFT.
    """
    N_sam = len(mu_components[0])
    autocorr_scipy = np.zeros(N_sam * 2 - 1)

    # 1 cm-1 = 2.99793e10 Hz  # https://wild.life.nctu.edu.tw/class/common/energy-unit-conv-table.html
    # 1 fs = 1e-15 s = 1e-15 1/Hz = 1e-15 / (1/2.99793e10 cm-1 ) = 1e-15 * 2.99793e10 cm = 2.99793e-5
    # dt_cm: "distance light travels in dt", so that fftfreq gives cm^-1
    # directly instead of Hz.
    dt_cm = dt_fs * 2.99793e-5

    for mu in mu_components:
        mu = np.asarray(mu, dtype=float)
        if use_gradient:
            mu = np.gradient(mu, dt_fs)  # d mu / d t
        mu = mu - np.mean(mu)  # remove average (important!)

        # ref https://stackoverflow.com/questions/643699/how-can-i-use-numpy-correlate-to-do-autocorrelation
        corr_result = correlate(mu, mu, mode='full', method='fft')
        # length 2*N_sam - 1, peaked at the center (zero lag);
        # goes from 0 to a large value back to 0.  .../\...
        autocorr_scipy += corr_result

    half = autocorr_scipy.size // 2  # floor division, index of zero lag
    cut = int(cut_autocorr / dt_fs)  # unit of cut: point
    cut = min(cut, half)  # guard against cut_autocorr > trajectory length
    autocorr_scipy = autocorr_scipy[(half - cut):(half + cut)]
    plot_debug(autocorr_scipy, verbose, 'before the window')

    if use_window:
        win = hann(2 * cut)
        autocorr_scipy = autocorr_scipy * win
        plot_debug(win, verbose, 'window')
        plot_debug(autocorr_scipy, verbose, 'after the window')

    N_sam_scipy = len(autocorr_scipy)
    # np.abs is needed because fft multiplies by exp(i omega t), a complex number.
    ir = np.abs(fft(autocorr_scipy)[0:N_sam_scipy // 2])
    freq_cm = fftfreq(N_sam_scipy, dt_cm)[:N_sam_scipy // 2]  # //: floor division

    if not use_gradient:
        ir = ir * freq_cm[:len(ir)] ** 2  # omega^2

    if normalize:
        ir = ir / np.max(ir)

    return freq_cm, ir


def get_ir_dipole():
    """This function obtains the IR intensity using the autocorrelation of dipole moments, or FFT of dipole
    moments."""
    import parameter as param
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
    use_gradient = True  # d mu / d t

    freq_range = getattr(param, 'freq_range', None)
    cut_autocorr = getattr(param, 'cut_autocorr', 1000)  # fs  # see MACE4IR paper

    # =========================
    # Load dipole data
    # =========================
    data = np.loadtxt(dipole_file)[:-1]
    # [:-1]: the length is a multiple of 2, probably FFT is slightly faster?

    mux = data[:, 1]
    muy = data[:, 2]
    muz = data[:, 3]

    ir_fig_verbose = getattr(param, 'ir_fig_verbose', False)
    freq_cm, ir = compute_ir_from_dipole_components(
        [mux, muy, muz], dt_fs,
        use_window=use_window, use_gradient=use_gradient,
        cut_autocorr=cut_autocorr, verbose=ir_fig_verbose,
    )

    # =========================
    # Save spectrum
    # =========================
    np.savetxt(output_data,
               np.column_stack((freq_cm, ir)),
               header="Frequency(cm^-1)  Intensity (AU)", fmt='%.6g')

    # =========================
    # Plot
    # =========================
    plt.figure(figsize=(7, 5))
    plt.plot(freq_cm, ir, linewidth=1.5)
    if freq_range is not None:
        plt.xlim(freq_range)
    plt.ylim([0, 1.05])
    plt.xlabel("Wavenumber (cm$^{-1}$)")
    plt.ylabel("Intensity")
    plt.tight_layout()
    plt.savefig(output_plot, dpi=300)
    if ir_fig_verbose:
        plt.show()


if __name__ == '__main__':
    get_ir_dipole()
