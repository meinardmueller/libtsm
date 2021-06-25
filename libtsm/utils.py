"""
Description: libtsm utility functions
Contributors: Sebastian Rosenzweig, Simon Schwär, Jonathan Driedger, Meinard Müller
License: The MIT license, https://opensource.org/licenses/MIT
This file is part of libtsm (https://www.audiolabs-erlangen.de/resources/MIR/2021-DAFX-AdaptivePitchShifting)
"""

import numpy as np
import scipy as sc
import scipy.signal
from typing import Tuple


def win(win_len, beta) -> np.ndarray:
    """
    Generates a sin^beta window.

    Parameters
    ----------
    win_len : int
             length of the window

    beta :   int
             Exponent of the window
    Returns
    -------
    w : np.ndarray [shape=(win_len, )]
        The desired window
    """
    w = np.sin((np.pi * np.arange(0, win_len)) / win_len) ** beta

    return w


def hps(x, ana_hop=256, win_length=1024, win_beta=2, Fs=22050, zero_pad=0, fil_len_harm=10, fil_len_perc=10,
        masking_mode='binary') -> Tuple[np.ndarray, np.ndarray]:
    """
    Harmonic - Percussive separation usign median filters.

    Parameters
    ----------
    x : np.ndarray [shape=(N, )], real - valued
        Signal to be transformed

    ana_hop : int
        hop size of the synthesis window

    win_length : int
        length of analysis and synthesis window for STFT

    win_beta : int
        exponent of sin^beta window

    Fs : int
        Sampling rate of the input audio signal x

    zero_pad : int
        For FFT. Number of zeros padded to the window to increase the fft size

    fil_len_harm: int
        Length of the median filter in time direction. A shorter filter makes it more likely
        that the signal is interpreted as harmonic.

    fil_len_perc: int
        Length of the median filter in frequency direction. A shorter filter makes it more likely
        that the signal is interpreted as percussive.

    masking_mode : either "binary" or "relative"
        Selects Harmonic Percussive separation masking mode (soft or binary masking)

    Returns
    -------
    x_harm : np.ndarray
        Harmonic Component of input signal x

    x_perc : np.ndarray
        Percussive Component of input signal x
    """
    # Pre calculations
    window = win(win_length, win_beta)

    if len(x.shape) == 1:
        x = x.reshape(-1, 1)

    num_of_chan = x.shape[1]

    #  harmonic-percussive separation
    x_harm = np.zeros(x.shape)  # Initialize output
    x_perc = np.zeros(x.shape)  # Initialize output

    for c in range(num_of_chan):  # loop over channels
        x_c = x[:, c]

        # stft
        spec, f, t = stft(x_c, ana_hop=ana_hop, win_length=win_length, win_beta=win_beta, Fs=Fs, num_of_frames=-1,
                          fft_shift=False, zero_pad=0)
        mag_spec = np.abs(spec)

        # harmonic-percussive separation
        mag_spec_perc = median_filter(mag_spec, fil_len_perc, 0)
        mag_spec_harm = median_filter(mag_spec, fil_len_harm, 1)

        if masking_mode == 'binary':
            mask_harm = mag_spec_harm > mag_spec_perc
            mask_perc = mag_spec_harm <= mag_spec_perc

        elif masking_mode == 'relative':
            mask_harm = mag_spec_harm / (mag_spec_harm + mag_spec_perc + np.finfo(float).eps)
            mask_perc = mag_spec_perc / (mag_spec_harm + mag_spec_perc + np.finfo(float).eps)

        else:
            raise Exception('masking mode must either be "binary" or "relative"!')

        spec_harm = mask_harm * spec
        spec_perc = mask_perc * spec

        # istft
        x_harm_c = istft(spec_harm, syn_hop=ana_hop, win_length=win_length, win_beta=win_beta, Fs=Fs, zero_pad=zero_pad,
                         num_of_iter=1, orig_sig_len=x.shape[0], restore_energy=False, fft_shift=False)
        x_perc_c = istft(spec_perc, syn_hop=ana_hop, win_length=win_length, win_beta=win_beta, Fs=Fs, zero_pad=zero_pad,
                         num_of_iter=1, orig_sig_len=x.shape[0], restore_energy=False, fft_shift=False)

        x_harm[:, c] = x_harm_c[:, 0]
        x_perc[:, c] = x_perc_c[:, 0]

    return x_harm, x_perc


def median_filter(X, filt_len, dim) -> np.ndarray:
    """
    Median filter implementation.

    Parameters
    ----------
    X : np.ndarray
        Spectrogram
    filt_len : int
        Median filter length
    dim : int
        Dimension in which median filter should be applied

    Returns
    -------
    Y : np.ndarray
        Median-filtered spectrogram

    """
    s = X.shape
    Y = np.zeros(s)

    if dim == 0:
        X_padded = np.vstack((np.zeros((int(np.floor(filt_len / 2)), s[1])),
                              X,
                              np.zeros((int(np.ceil(filt_len / 2)), s[1]))))
        for i in range(s[0]):
            Y[i, :] = np.median(X_padded[i:i + filt_len, :], axis=0)

    elif dim == 1:
        X_padded = np.hstack((np.zeros((s[0], int(np.floor(filt_len / 2)))),
                              X,
                              np.zeros((s[0], int(np.ceil(filt_len / 2))))))
        for i in range(s[1]):
            Y[:, i] = np.median(X_padded[:, i:i + filt_len], axis=1)

    else:
        raise Exception("Invalid div!")

    return Y


def stft(x, ana_hop=2048, win_length=4096, win_beta=2, Fs=22050, num_of_frames=-1, fft_shift=False, zero_pad=0) -> \
        Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Computes the Short-Time Fourier Transform (STFT) of the input audio signal.

    Parameters
    ----------
    x : np.ndarray, real-valued
        Signal to be transformed

    ana_hop : int or np.ndarray
        hop size of the analysis window

    win_length : int
        length of analysis window for STFT

    win_beta : int
        exponent of sin^beta window

    Fs : int
        Sampling rate of the input audio signal x

    num_of_frames : int
        Fixes the number of FFT frames to be computed

    fft_shift: bool
        For FFT. When True, applies a circular shift to each frame of half its length, prior to computing the FFT

    zero_pad : int
        For FFT. Number of zeros padded to the window to increase the fft size

    Returns
    -------
    X : np.ndarray [shape=(K, M + 1)], complex-valued
        The discrete short-time Fourier transform

    f : np.ndarray [shape=(K, )], real-valued
        Center frequencies of all Fourier bins given in Hertz

    t : np.ndarray [shape=(M+1, )], real-valued
        Time instances where the respective Fourier spectra were computed, given in seconds

    """
    # Pre-calculations
    window = win(win_length, win_beta)

    # Zero-pad the window
    w = np.concatenate((np.zeros(int(np.floor(zero_pad / 2))), window, np.zeros(int(np.floor(zero_pad / 2)))))
    win_len = int(len(w))
    win_len_half = np.around(win_len / 2).astype(int)

    max_ana_hop = int(np.max(ana_hop))

    if len(x.shape) == 1:
        x = x.reshape(-1, 1)

    # Pad the audio to center the windows and to avoid problems at the end
    x_padded = np.vstack((np.zeros((win_len_half, 1)), x, np.zeros((win_len+max_ana_hop+1, 1))))

    # In case ana_hop is a scalar, sample the window positions evenly in the input signal
    if np.isscalar(ana_hop):
        if num_of_frames < 0:
            num_of_frames = int(np.floor((len(x_padded) - win_len)/ana_hop + 1))

        win_pos = np.arange(num_of_frames).astype(int) * ana_hop
    else:
        if num_of_frames < 0:
            num_of_frames = len(ana_hop)

        win_pos = ana_hop[0:num_of_frames].astype(int)

    # Spectrogram calculation
    spec = np.zeros((win_len_half + 1, num_of_frames), dtype=complex)

    for i in range(num_of_frames):
        xi = x_padded[win_pos[i]:win_pos[i] + win_len] * w.reshape(-1, 1)

        if fft_shift == 1:
            xi = np.fft.fftshift(xi)

        Xi = np.fft.fft(xi, axis=0)

        spec[:, i] = Xi[0:win_len_half + 1, 0]

    # Axis calculation
    t = win_pos / Fs
    f = np.arange(0, win_len_half + 1) * Fs / win_len

    return spec, f, t


def istft(spec, syn_hop=2048, win_length=4096, win_beta=2, Fs=22050, zero_pad=0, num_of_iter=1, orig_sig_len=-1,
          restore_energy=False, fft_shift=False) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Computes the 'inverse' Short Time Fourier Transform, according to the paper "Signal Estimation from Modified
    Short-Time Fourier Transform" by Griffin and Lim.

    Parameters
    ----------
    spec : np.ndarray [shape=(K,M+1)] , complex-valued
        A complex spectrogram generated by STFT.

    syn_hop : int
        hop size of the synthesis window

    win_length : int
        length of synthesis window for ISTFT

    win_beta : int
        exponent of sin^beta window

    Fs : int
        sampling rate

    zero_pad : int
        For IFFT. Number of zeros padded to the window to increase the fft size

    num_of_iter : int
        number of iterations for synthesis

    orig_sig_len : int
        Original length of the audio signal such that the output can be trimmed accordingly, in samples

    restore_energy : bool
        For IFFT. When True, rescales every windowed synthesis frame to compensate for synthesis energy leakage

    fft_shift : bool
        For IFFT. When True, applies a circular shift to each frame of half its length, prior to computing the FFT

    Returns
    -------
    y : np.ndarray [shape=(L,1)], real - valued
        The time-domain signal.
    """

    # Pre-calculations
    num_of_frames = spec.shape[1]

    # First iteration
    Y_i = spec
    y_i = lsee_mstft(Y_i, syn_hop=syn_hop, win_length=win_length, win_beta=win_beta, zero_pad=zero_pad,
                     restore_energy=restore_energy, fft_shift=fft_shift)

    # Remaining iterations
    for j in range(1, num_of_iter):
        Y_i = np.abs(spec) * np.exp(1j*np.angle(stft(y_i, ana_hop=syn_hop, win_length=win_length, win_beta=win_beta,
                                                     Fs=Fs, num_of_frames=num_of_frames, fft_shift=fft_shift,
                                                     zero_pad=zero_pad)[0]))
        y_i = lsee_mstft(Y_i, syn_hop=syn_hop, win_length=win_length, win_beta=win_beta, zero_pad=zero_pad,
                         restore_energy=restore_energy, fft_shift=fft_shift)

    y = y_i

    # If the original Length of the signal is known, also remove the zero padding at the end
    if orig_sig_len > 0:
        y = y[:orig_sig_len]

    return y


def lsee_mstft(X, syn_hop=2048, win_length=4096, win_beta=2, zero_pad=0, restore_energy=0, fft_shift=0) -> np.ndarray:
    """
    Computes the 'inverse' Short Time Fourier Transform (ISTFT) using the Griffin Lim procedure.

    Parameters
    ----------
    X : np.ndarray [shape=(K,M+1)] , complex-valued
        A complex spectrogram generated by STFT.

    syn_hop : int
        hop size of the synthesis window

    win_length : int
        length of analysis and synthesis window for STFT

    win_beta : int
        exponent of sin^beta window

    zero_pad : int
        For IFFT. Number of zeros padded to the window to increase the fft size

    restore_energy : bool
        For IFFT. When True, rescales every windowed synthesis frame to compensate for synthesis energy leakage

    fft_shift : bool
        For IFFT. When True, applies a circular shift to each frame of half its length, prior to computing the FFT

    Returns
    -------
    x: np.ndarray [shape=(L,1)], real-valued
        The time-domain signal.
    """

    # Pre-calculations
    window = win(win_length, win_beta)
    w = np.concatenate((np.zeros(int(np.floor(zero_pad / 2))), window, np.zeros(int(np.floor(zero_pad / 2)))))
    win_len = int(len(w))
    win_len_half = np.around(win_len / 2).astype(int)
    num_of_frames = X.shape[1]
    win_pos = np.arange(0, num_of_frames).astype(int) * syn_hop
    signal_length = win_pos[-1] + win_len

    x = np.zeros((signal_length, 1))  # re-synthesized signal
    ow = np.zeros((signal_length, 1))  # sum of overlapping windows

    for i in range(num_of_frames):
        curr_spec = X[:, i]

        # add the conjugate complex symmetric upper half of the spectrum
        Xi = np.concatenate((curr_spec, np.conj(curr_spec[-2:0:-1])))
        xi = np.real(np.fft.ifft(Xi, axis=0))

        if fft_shift == 1:
            xi = np.fft.fftshift(xi)

        xiw = xi * w

        if restore_energy == 1:
            xi_energy = np.sum(np.abs(xi))
            xi_w_energy = np.sum(np.abs(xiw))
            xiw = xiw * (xi_energy/(xi_w_energy+np.finfo(float).eps))

        x[win_pos[i]:win_pos[i] + win_len, 0] += xiw
        ow[win_pos[i]:win_pos[i] + win_len, 0] += w**2

    ow[ow < 10**-3] = 1  # avoid potential division by zero
    x = x / ow

    # knowing the zeropads that were added in the stft computation, we can remove them again now. But since we do not
    # know exactly how many zeros were padded at the end of the signal, it is only safe to remove winLenHalf zeros.
    x = x[win_len_half:-win_len_half, :]

    return x


def find_peaks(X) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Finds peaks on spectrum X. An index in X is considered a peak if its value is the largest among its four nearest
    neighbours.

    Parameters
    ----------
    X : np.ndarray [shape=(K, )] , complex-valued
        An FFT vector.

    Returns
    -------
    peaks : np.ndarray [shape=(P, )] , real-valued
        Vector with P peaks found
    """
    mag_spec = np.abs(X)
    mag_spec_padded = np.hstack((np.zeros(2), mag_spec, np.zeros(2)))
    peaks = np.where((mag_spec_padded[4:] < mag_spec_padded[2:-2]) &
                     (mag_spec_padded[3:-1] < mag_spec_padded[2:-2]) &
                     (mag_spec_padded[1:-3] < mag_spec_padded[2:-2]) &
                     (mag_spec_padded[0:-4] < mag_spec_padded[2:-2]))[0]

    infl_region_start = np.zeros(peaks.shape, dtype=int)
    infl_region_end = np.zeros(peaks.shape, dtype=int)

    if peaks.size == 0:
        return peaks, infl_region_start, infl_region_end

    infl_region_start[0] = 0
    infl_region_start[1:] = np.ceil((peaks[1:] + peaks[0:-1])/2)
    infl_region_end[0:-1] = infl_region_start[1:]
    infl_region_end[-1] = len(infl_region_end)

    return peaks, infl_region_start, infl_region_end


def cross_corr(x, y, win_len) -> np.ndarray:
    """
    Computes cross correlation between signals x and y over a window of size win_len.

    Parameters
    ----------
    x : np.ndarray [shape=(N, )], real or complex - valued
        Signal to be cross-correlated

    y : np.ndarray [shape=(N, )], real or complex - valued
        Signal to be cross-correlated

    win_len : int
        Cross correlation window, in samples

    Returns
    -------
    y : np.ndarray [shape=(2N+1,)], real - valued
        Crosscorrelated signal
    """
    # cross correlation is essentially the same as convolution with the first signal being reverted. In principle, we
    # also need to take the complex conjugate of the reversed x, but since audio signals are real valued, we can skip
    # this operation.
    cc = np.convolve(np.flip(x), y)

    # restrict the cross correlation result to just the relevant values
    # Values outside of this range are related to deltas bigger or smaller than our tolerance values.
    cc = cc[win_len-1:-(win_len-1)]

    return cc


def normalize_length(xf, length) -> np.ndarray:
    """
    Adjusts the length of signal xf to variable "length".

    Parameters
    ----------
    xf : np.ndarray [shape=(N, )], real or complex - valued
        Signal to be processed

    length : int
        Signal to be cross-correlated

    Returns
    -------
    y : np.ndarray [shape=(length, )], real or complex - valued
        Signal with modified length
    """

    if len(xf[:, 0]) < length:
        pad_len = length - len(xf[:, 0])
        y = np.concatenate((xf, np.zeros((pad_len, xf.shape[1]))), axis=0)
    else:
        y = xf[0:length, :]
    return y


def modify_spectral_envelope(x, y, ana_hop=64, win_length=1024, win_beta=1, Fs=22050,  filt_len=24) -> np.ndarray:
    """
    Complement to the pitch shifting algorithm, that modifies the formants of the
    pitch-shifted signal to match them with those of the original signal.

    Parameters
    ----------
    x : np.ndarray [shape=(N, )], real - valued
        Original input signal

    y : np.ndarray [shape=(N, )], real - valued
        Pitch-shifted signal

    ana_hop : int
        hop size of the STFT analysis and synthesis window

    win_length : int
        length of the analysis and synthesis window for STFT

    win_beta : int
        exponent of sin^beta window

    Fs : int
        Sampling rate of audio signals x and y

    filt_len : int
        number of samples for envelope modifying function

    Returns
    -------
    y_spec_env_X : np.ndarray [shape=(N,)], real - valued
        Pitch-shifted signal with modified spectral envelope
    """

    if len(x.shape) == 1:
        x = x.reshape(-1, 1)

    if len(y.shape) == 1:
        y = y.reshape(-1, 1)

    num_of_chan = x.shape[1]
    y_spec_env_x = np.zeros(y.shape)

    for c in range(num_of_chan):
        x_c = x[:, c]
        y_c = y[:, c]

        # stft
        X, _, _ = stft(x_c, ana_hop=ana_hop, win_length=win_length, win_beta=win_beta, Fs=Fs, num_of_frames=-1,
                       fft_shift=False, zero_pad=0)
        Y, _, _ = stft(y_c, ana_hop=ana_hop, win_length=win_length, win_beta=win_beta, Fs=Fs, num_of_frames=-1,
                       fft_shift=False, zero_pad=0)

        # Compute spectral envelopes
        env_X = comp_env(X, filt_len)
        env_Y = comp_env(Y, filt_len)
        Y_spec_env_X = np.multiply(np.divide(Y, env_Y), env_X)

        # istft
        y_spec_env_x[:, c] = istft(Y_spec_env_X, syn_hop=ana_hop, win_length=win_length, win_beta=win_beta, Fs=Fs,
                                   zero_pad=0, num_of_iter=1, orig_sig_len=len(x), restore_energy=False,
                                   fft_shift=False)[:, 0]

    return y_spec_env_x


def comp_env(X, filt_len) -> np.ndarray:
    """
    Computes the envelope of a given signal spectrum.

    Parameters
    ----------
    X : np.ndarray [shape=(K,M+1)] , complex-valued
        A complex spectrogram

    filt_len : int
        Length of the convolution window

    Returns
    -------
    env : np.ndarray [shape=(K,M+1)] , real-valued
        Spectral Envelope
    """
    kern = win(filt_len, 2)  # Hann Window
    kern.shape = (-1, 1)  # Turn the window into a 2d array
    env = sc.signal.convolve2d(np.abs(X), kern, mode='same')  # not equivalent to Matlab's conv2()
    env = np.divide(env, np.finfo(float).eps + np.tile(np.max(env, axis=0), (env.shape[0], 1)))  # Normalization
    env = np.multiply(env, np.tile(np.max(np.abs(X), axis=0), (np.abs(X).shape[0], 1)))  # Scaling
    env[env < 0.01] = 0.01

    return env
