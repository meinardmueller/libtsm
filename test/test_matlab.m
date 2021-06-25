%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test script for executing Matlab functions
% Contributors: Sebastian Rosenzweig, Simon Schwaer, Jonathan Driedger, Meinard Mueller
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
