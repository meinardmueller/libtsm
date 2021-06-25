"""
Description: Tests for numerical comparison of Matlab implementation and libtsm
Contributors: Sebastian Rosenzweig, Simon Schwär, Jonathan Driedger, Meinard Müller
License: The MIT license, https://opensource.org/licenses/MIT
This file is part of libtsm (https://www.audiolabs-erlangen.de/resources/MIR/2021-DAFX-AdaptivePitchShifting)
"""

import numpy as np
import scipy.io
import os.path
import soundfile as sf
import libtsm


# this file is needed for comparing the Matlab outputs to the libtsm outputs
matlab_output = "output/matlab.mat"
assert os.path.exists(matlab_output), "Please execute test_matlab.m first!"
ml = scipy.io.loadmat(matlab_output)  # load matlab output

audio_file = "data/CastanetsViolin.wav"
x, Fs = sf.read(audio_file)
assert Fs == 22050


def test_x(skip_numerical):
    assert ml["x"][:, 0].shape == x.shape


# 1. Utils #############################################################################################################
def test_win(skip_numerical):
    window = libtsm.win(1024, 2)

    assert ml["window"][:, 0].shape == window.shape
    assert np.allclose(ml["window"][:, 0], window)


def test_stft(skip_numerical):
    Y, f, t = libtsm.stft(x, ana_hop=2048, win_length=4096, win_beta=2)

    assert ml["Y"].shape == Y.shape
    assert ml["f"][0, :].shape == f.shape
    assert ml["t"][0, :].shape == t.shape
    assert np.allclose(ml["Y"], Y)
    assert np.allclose(ml["f"][0, :], f)
    assert np.allclose(ml["t"][0, :], t)


def test_istft(skip_numerical):
    x_i = libtsm.istft(ml["Y"], syn_hop=2048, win_length=4096, win_beta=2)

    assert ml["xI"].shape == x_i.shape
    assert np.allclose(ml["xI"], x_i)


def test_modify_spectral_envelope(skip_numerical):
    y_spec_env = libtsm.modify_spectral_envelope(x, x)

    assert ml["ySpecEnv"].shape == y_spec_env.shape
    assert np.allclose(ml["ySpecEnv"], y_spec_env)


def test_hps(skip_numerical):
    x_harm, x_perc = libtsm.hps(x, Fs=Fs)

    assert ml["xHarm"].shape == x_harm.shape
    assert ml["xPerc"].shape == x_perc.shape
    assert np.allclose(ml["xHarm"], x_harm)
    assert np.allclose(ml["xPerc"], x_perc)


# 2. TSM ###############################################################################################################
def test_ola_tsm(skip_numerical):
    alpha = 1.8  # time-stretch factor
    y_ola = libtsm.wsola_tsm(x, alpha, syn_hop=128, win_length=256, win_beta=2, tol=0)

    assert ml["yOLA"].shape == y_ola.shape
    assert np.allclose(ml["yOLA"], y_ola)


def test_wsola_tsm(skip_numerical):
    alpha = 1.8  # time-stretch factor
    y_wsola = libtsm.wsola_tsm(x, alpha)

    assert ml["yWSOLA"].shape == y_wsola.shape
    assert np.allclose(ml["yWSOLA"], y_wsola)


def test_pv_tsm(skip_numerical):
    y_pv = libtsm.pv_tsm(x, 1.8)

    assert ml["yPV"].shape == y_pv.shape
    assert np.allclose(ml["yPV"], y_pv, atol=10**(-2), rtol=10**(-2))  # slight deviations due to rounding errors


def test_pvpl_tsm(skip_numerical):
    y_pvpl = libtsm.pv_tsm(x, 1.8, phase_locking=True)

    assert ml["yPVpl"].shape == y_pvpl.shape
    assert np.allclose(ml["yPVpl"], y_pvpl)


def test_hps_tsm(skip_numerical):
    alpha = 1.8
    y_hp = libtsm.hps_tsm(x, alpha=alpha)

    assert ml["yHP"].shape == y_hp.shape
    assert np.allclose(ml["yHP"], y_hp)


def test_pv_int_tsm(skip_numerical):
    y_pv_int = libtsm.pv_int_tsm(x, alpha=4)

    assert ml["yPvInt"].shape == y_pv_int.shape
    assert np.allclose(ml["yPvInt"], y_pv_int)


def test_two_step_tsm(skip_numerical):
    alpha = 1.8
    y_two_step = libtsm.two_step_tsm(x, alpha)

    assert ml["yTwoStep"].shape == y_two_step.shape
    assert np.allclose(ml["yTwoStep"], y_two_step)


def test_nonlinear_tsm(skip_numerical):
    audio_file1 = "data/BeethovenOrchestra.wav"
    x1, Fs1 = sf.read(audio_file1)
    assert Fs1 == 22050

    # extract anchor points from .MAT file
    mat_file = "data/BeethovenAnchorpoints.mat"
    mat = scipy.io.loadmat(mat_file)
    anchor_points = mat["anchorpoints"] - 1  # substract 1 for Python version

    # HPS-TSM using anchorpoints to synchronize orchestra with Piano file
    y_sync = libtsm.hps_tsm(x1, anchor_points)

    assert ml["ySync"].shape == y_sync.shape
    assert np.allclose(ml["ySync"], y_sync, atol=10**(0), rtol=10**(0))  # partial deviations due to rounding errors
