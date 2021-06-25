%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: demoTSMtoolbox
% Date: 12-2013
% Programmer: Jonathan Driedger
% http://www.audiolabs-erlangen.de/resources/MIR/TSMtoolbox/
%
% This is the demo script which illustrates the main functionalities of the
% 'TSM toolbox'. For a detailed description we refer to  [DM14] (see
% References below).
%
% The script proceeds in the following steps:
%   1. It loads a wav file.
%   2. It applies different TSM algorithms to the loaded signal with a
%      fixed constant stretching factor.
%   3. It visualizes the TSM results.
%   4. It writes the computed TSM results as wav files to the hard drive.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reference:
% If you use the 'TSM toolbox' please refer to:
% [DM14] Jonathan Driedger, Meinard Mueller
%        TSM Toolbox: MATLAB Implementations of Time-Scale Modification
%        Algorithms
%        Proceedings of the 17th International Conference on Digital Audio
%        Effects, Erlangen, Germany, 2014.
%
% License:
% This file is part of 'TSM toolbox'.
%
% MIT License
% 
% Copyright (c) 2021 Jonathan Driedger, Meinard Mueller, International Audio
% Laboratories Erlangen, Germany.
% 
% We thank the German Research Foundation (DFG) for various research grants 
% that allow us for conducting fundamental research in music processing.
% The International Audio Laboratories Erlangen are a joint institution of 
% the Friedrich-Alexander-Universitaet Erlangen-Nuernberg (FAU) and
% Fraunhofer Institute for Integrated Circuits IIS.
% 
% Permission is hereby granted, free of charge, to any person obtaining a 
% copy of this software and associated documentation files (the 
% "Software"), to deal in the Software without restriction, including 
% without limitation the rights to use, copy, modify, merge, publish, 
% distribute, sublicense, and/or sell copies of the Software, and to permit 
% persons to whom the Software is furnished to do so, subject to the 
% following conditions:
% 
% The above copyright notice and this permission notice shall be included 
% in all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
% IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
% CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
% TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
% SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('../matlab/'))
pathData = '../data/';
pathOutput = '../output/';

filename = 'CastanetsViolin.wav';
% filename = 'Bongo.wav';
% filename = 'DrumSolo.wav';
% filename = 'Glockenspiel.wav';
% filename = 'Jazz.wav';
% filename = 'Pop.wav';
% filename = 'SingingVoice.wav';
% filename = 'Stepdad.wav';
% filename = 'SynthMono.wav';
% filename = 'SynthPoly.wav';

alpha = 1.8; % constant stretching factor

warning('OFF','MATLAB:audiovideo:audiowrite:dataClipped');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Load the audio signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x,fsAudio] = audioread([pathData filename]);
%x = sin(2*pi*440*[0:1/fsAudio:(length(x)/fsAudio)-(1/fsAudio)])';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Auxiliary functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% window function
window = win(1024, 2);

% STFT
[Y, f, t] = stft(x);

% ISTFT
xI = istft(Y);

% Harmonic Percussive Separation
parHps = [];
[xHarm, xPerc] = hpSep(x, parHps);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. TSM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OLA
paramOLA.tolerance = 0;
paramOLA.synHop = 128;
paramOLA.win = win(256,2); % hann window
yOLA = wsolaTSM(x,alpha,paramOLA);

% WSOLA
yWSOLA = wsolaTSM(x,alpha);

% Phase Vocoder
yPV = pvTSM(x,alpha);

% Phase Vocoder with identity phase locking
paramPVpl.phaseLocking = 1;
yPVpl = pvTSM(x,alpha,paramPVpl);

% TSM based on HPSS
yHP = hpTSM(x,alpha);

% efficient TSM implementation
yPvInt = pvIntTSM(x, 4);

% Two-step TSM
yTwoStep = twoStepTSM(x, alpha);

% non-linear TSM
filenameWav2 = 'BeethovenOrchestra.wav';
filenameAP = 'BeethovenAnchorpoints.mat';
[x2,fsAudio] = audioread([pathData filenameWav2]);
mat = load([pathData filenameAP]);
anchorpoints = mat.anchorpoints;

ySync = hpTSM(x2,anchorpoints);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Pitch-Shifting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pitch-shift TSM
n = 100;
yPitchShift = pitchShiftViaTSM(x, n);

% Modify spectral envelope
ySpecEnv = modifySpectralEnvelope(x, x);


% save all variables
save([pathOutput 'matlab.mat'], 'x', 'window', 'Y', 'f', 't', 'xI', 'yOLA', 'yWSOLA', 'yPV', 'yPVpl', 'yHP', 'xHarm', 'xPerc', 'yPvInt', 'yTwoStep', 'yPitchShift', 'ySpecEnv', 'ySync');


