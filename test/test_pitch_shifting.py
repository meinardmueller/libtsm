"""
Description: tests of libtsm pitch-shifting functions
Contributors: Sebastian Rosenzweig, Simon Schwär, Jonathan Driedger, Meinard Müller
License: The MIT license, https://opensource.org/licenses/MIT
This file is part of libtsm (https://www.audiolabs-erlangen.de/resources/MIR/2021-DAFX-AdaptivePitchShifting)
"""

import numpy as np
import libtsm
import librosa


def test_fixed_pitch_shift_even():
    fn_in = 'data/three_sinusoidals.wav'
    x, Fs = librosa.load(fn_in, sr=22050)

    p = 1200  # cents
    y_psf1 = libtsm.pitch_shift(x, p, order='res-tsm')[:, 0]
    y_psf2 = libtsm.pitch_shift(x, p, order='tsm-res')[:, 0]
    y_psfo = libtsm.pitch_shift_original(x, p)[:, 0]

    assert len(x) == len(y_psf1)
    assert len(x) == len(y_psf2)
    assert len(x) == len(y_psfo)


def test_fixed_pitch_shift_odd():
    fn_in = 'data/three_sinusoidals.wav'
    x, Fs = librosa.load(fn_in, sr=22050)

    p = 1231  # cents
    y_psf1 = libtsm.pitch_shift(x, p, order='res-tsm')[:, 0]
    y_psf2 = libtsm.pitch_shift(x, p, order='tsm-res')[:, 0]
    y_psfo = libtsm.pitch_shift_original(x, p)[:, 0]

    assert len(x) == len(y_psf1)
    assert len(x) == len(y_psf2)
    assert len(x) == len(y_psfo)


def test_adaptive_pitch_shift1():
    # input signal
    sig_len = 10  # seconds
    Fs = 22050  # Hz
    t_sine = np.arange(0, sig_len, 1 / Fs)
    F_sine = 440
    sine = np.sin(2 * np.pi * F_sine * t_sine)

    # time-varying pitch-shift (sinusoidal)
    F_shift = 1  # Hz
    p = np.sin(2 * np.pi * F_shift * t_sine) * 200

    # pitch-shifting
    siren1 = libtsm.pitch_shift(sine, p, t_p=t_sine, order='res-tsm')[:, 0]
    siren2 = libtsm.pitch_shift(sine, p, t_p=t_sine, order='tsm-res')[:, 0]

    assert len(sine) == len(siren1)
    assert len(sine) == len(siren2)


def test_adaptive_pitch_shift2():
    # Load Audio File
    fn_in = 'data/three_sinusoidals.wav'
    x, Fs = librosa.load(fn_in, sr=22050)

    # Adaptive Pitch-Shifting (Figure 3)
    t = np.arange(0, len(x) / Fs, 1 / Fs)  # sec
    N = len(t)
    t_1 = t[0:N // 3]
    t_2 = t[N // 3:2 * N // 3]
    t_3 = t[2 * N // 3:]

    p = np.concatenate((np.zeros(len(t_1)), 800 * np.sin(2 * np.pi * 1 * t_2), np.linspace(0, 1200, len(t_3))))  # cents

    y_psa1 = libtsm.pitch_shift(x, p, t, order='res-tsm')[:, 0]
    y_psa2 = libtsm.pitch_shift(x, p, t, order='tsm-res')[:, 0]

    assert len(x) == len(y_psa1)
    assert len(x) == len(y_psa2)


def test_adaptive_pitch_shift3():
    Fs = 48000
    t = np.linspace(0, 5, int(5*Fs))
    x = 0.5 * np.sin(2 * np.pi * 220 * t)

    p = np.array([0, 100, 50, -100, 0])
    t_p = np.array([0, 1, 2, 3, 4])

    y1 = libtsm.pitch_shift(x, p, t_p, Fs=Fs, order='res-tsm')[:, 0]
    y2 = libtsm.pitch_shift(x, p, t_p, Fs=Fs, order='tsm-res')[:, 0]

    assert len(x) == len(y1)
    assert len(x) == len(y2)
