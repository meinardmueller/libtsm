"""
Description: libtsm pitch-shifting functions
Contributors: Sebastian Rosenzweig, Simon Schwär, Jonathan Driedger, Meinard Müller
License: The MIT license, https://opensource.org/licenses/MIT
This file is part of libtsm (https://www.audiolabs-erlangen.de/resources/MIR/2021-DAFX-AdaptivePitchShifting)
"""

import numpy as np
import scipy as sc
import scipy.signal
import scipy.interpolate
from fractions import Fraction as frac
from .tsm import hps_tsm
from .utils import normalize_length


def pitch_shift_original(x, n, Fs=22050) -> np.ndarray:
    """
    Pitch modification algorithm via time-scale modification. The input signal is stretched via TSM and then resampled.
    The code closely follows the Matlab implementation.

    Parameters
    ----------
    x : np.ndarray [shape=(N, )], real - valued
        Signal to be transformed

    n : int
        Amount of pitch shifting to be applied, given in cents. Positive n indicates pitch rising, negative n a pitch
        lowering

    Fs : int
        Sampling rate of the input audio signal x

    Returns
    -------
    y : np.ndarray [shape=(L,1)], real - valued
        The time-scale modified output signal
    """

    if len(x.shape) == 1:
        x = x.reshape(-1, 1)

    alpha = np.power(np.power(2, 1 / 12), (n / 100))
    y_tsm = hps_tsm(x, alpha, Fs=Fs)
    const = Fs / np.around(alpha * Fs)
    f = frac(str(const)).limit_denominator(100)
    p = f.numerator
    q = f.denominator
    x_f = sc.signal.resample_poly(y_tsm, int(p), int(q), axis=0)  # deviations from Matlab's resample()
    y = normalize_length(x_f, x.shape[0])

    return y


def pitch_shift(x, p, t_p=None, Fs=22050, order="res-tsm", **kwargs) -> np.ndarray:
    """
    (Non-linear) pitch-shifting via time-scale modification and resampling.

    Parameters
    ----------
    x : np.ndarray [shape=(N, )], real - valued
        Signal to be transformed

    p : float or np.ndarray [shape=(M,)], real - valued
        Amount of pitch shifting to be applied, given in cents. Positive p indicates pitch rising, negative p a pitch
        lowering.

    t_p : np.ndarray [shape=(M,)], real - valued
        Array of time instances in seconds for adaptive pitch shifting, same length as p. If t==None, a fixed
        pitch-shift is assumed.

    Fs : int
        Sampling rate of the input audio signal x

    order : Order of TSM and resampling, either "res-tsm" or "tsm-res".

    **kwargs : Parameters for hps_tsm

    Returns
    -------
    y : np.ndarray [shape=(L,1)], real - valued
        The time-scale modified output signal
    """

    if len(x.shape) == 1:
        x = x.reshape(-1, 1)

    if not np.isscalar(p):
        if t_p is None:
            raise Exception("t must be specified if p is an array!")
        if len(p) != len(t_p):
            raise Exception("t must have the same length as p!")

    t_x = np.linspace(0, (len(x) - 1) / Fs, len(x))

    # account for sign change when order of resampling and TSM is exchanged
    if order == "res-tsm":
        alpha = 2 ** (-p / 1200)
    elif order == "tsm-res":
        alpha = 2 ** (p / 1200)
    else:
        raise Exception("Order must be either res-tsm or tsm-res!")

    # convert pitch shift in cents to (non-linear) time-stretch function tau
    if np.isscalar(p):
        tau = np.array([[0, 0], [x.shape[0] - 1, x.shape[0] * alpha - 1]]) / Fs  # given in seconds
    else:
        # compute tau
        tau = np.zeros((len(alpha), 2))
        tau[:, 0] = t_p

        for i in range(1, len(alpha)):
            tau[i, 1] = alpha[i] / Fs + tau[i - 1, 1]

    # Pitch-shifting
    if order == "res-tsm":
        # (Non-linear) Resampling
        fi = sc.interpolate.interp1d(tau[:, 0], tau[:, 1], kind='linear', fill_value="extrapolate")
        time_input = fi(t_x)
        fi = sc.interpolate.interp1d(time_input, x[:, 0], kind='cubic', fill_value="extrapolate")
        t_res = np.arange(0, tau[-1, 1] + 1 / Fs, 1 / Fs)
        y_ps = fi(t_res)

        tau_inv = np.hstack((time_input.reshape(-1, 1), t_x.reshape(-1, 1)))
        anchor_points = np.round(tau_inv * Fs).astype(int)
        anchor_points = anchor_points[np.unique(anchor_points[:, 0],
                                                return_index=True)[1], :]  # only keep unique indices

        # Time-Scale Modification
        y_ps = hps_tsm(y_ps, anchor_points, Fs=Fs, **kwargs)

    elif order == "tsm-res":
        # compute anchor points
        anchor_points = np.round(tau * Fs).astype(int)
        anchor_points = anchor_points[np.unique(anchor_points[:, 1],
                                                return_index=True)[1], :]  # only keep unique indices

        # Time-Scale Modification
        y_tsm = hps_tsm(x, anchor_points, Fs=Fs, **kwargs)

        # (Non-linear) resampling
        time_output = np.linspace(0, (y_tsm.shape[0] - 1) / Fs, y_tsm.shape[0])
        fi = sc.interpolate.interp1d(tau[:, 1], tau[:, 0], kind='linear', fill_value="extrapolate")
        time_input = fi(time_output)
        fi = sc.interpolate.interp1d(time_input, y_tsm[:, 0], kind='cubic', fill_value="extrapolate")
        y_ps = fi(t_x)

    return y_ps.reshape(-1, 1)

