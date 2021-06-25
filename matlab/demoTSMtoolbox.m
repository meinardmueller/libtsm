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
% 1. load the audio signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x,fsAudio] = audioread([pathData filename]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. compute TSM results
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

% elastique
% To execute elastique, you will need an access id from
% http://www.sonicapi.com. Furthermore, you need to download 'curl'
% from http://curl.haxx.se/download.html. For further information, see the
% function header of the file elastiqueTSM.m.

% paramElast.fsAudio = fsAudio;
% yELAST = elastiqueTSM(x,alpha,paramElast);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. visualize TSM results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Original
timeRange =  [0.0, 0.5];
paramVis.timeRange = timeRange;
paramVis.title = 'Original';
visualizeWav(x,paramVis);

timeRange = timeRange * alpha;
paramVis.timeRange = timeRange;

% OLA
paramVis.title = 'OLA';
visualizeWav(yOLA,paramVis);

% WSOLA
paramVis.title = 'WSOLA';
visualizeWav(yWSOLA,paramVis);

% Phase Vocoder
paramVis.title = 'Phase Vocoder';
visualizeWav(yPV,paramVis);

% Phase Vocoder with identity phase locking
paramVis.title = 'Phase Vocoder with identity phase locking';
visualizeWav(yPVpl,paramVis);

% TSM based on HPSS
paramVis.title = 'HP-TSM';
visualizeWav(yHP,paramVis);

% elastique
% paramVis.title = 'elastique';
% visualizeWav(yELAST,paramVis);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. write TSM results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OLA
audiowrite([pathOutput filename(1:end-4) '_' sprintf('%0.2f',alpha) '_OLA.wav'],...
            yOLA,fsAudio);

% WSOLA
audiowrite([pathOutput filename(1:end-4) '_' sprintf('%0.2f',alpha) '_WSOLA.wav'],...
            yWSOLA,fsAudio);

% Phase Vocoder
audiowrite([pathOutput filename(1:end-4) '_' sprintf('%0.2f',alpha) '_PV.wav'],...
            yPV,fsAudio);

% Phase Vocoder with identity phase locking
audiowrite([pathOutput filename(1:end-4) '_' sprintf('%0.2f',alpha) '_PVpl.wav'],...
            yPVpl,fsAudio);

% TSM based on HPSS
audiowrite([pathOutput filename(1:end-4) '_' sprintf('%0.2f',alpha) '_HP.wav'],...
            yHP,fsAudio);

% elastique
% audiowrite([pathOutput filename(1:end-4) '_' sprintf('%0.2f',alpha) '_ELAST.wav'],...
%             yELAST,fsAudio);
