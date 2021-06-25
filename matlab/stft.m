function [spec,f,t,sideinfo] = stft(x,parameter,sideinfo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: stft
% Date: 12-2013
% Programmer: Jonathan Driedger
% http://www.audiolabs-erlangen.de/resources/MIR/TSMtoolbox/
%
% Computes the short-time Fourier transform (stft) of the input audio
% signal.
%
% Input:    x               single-channel signal.
%           parameter.
%               anaHop      either the constant hop size of the analysis
%                           window or a vector of analysis positions in the
%                           input signal.
%               win         the analysis window used for windowing the
%                           input signal.
%               zeroPad     number of zeros that should be padded to the
%                           window to increase the fft size and therefore
%                           the frequency resolution.
%               fsAudio     the sampling rate of the input audio signal x.
%               numOfFrames can be set to fix the number of spectra that
%                           should be computed.
%               fftShift    can be set to 1 to apply a circular shift of
%                           samples to each frame by half its length prior
%                           to the application of the fft.
%
% Output:   spec            complex spectrogram
%           f               vector of the center frequencies of all Fourier
%                           bins given in Hertz.
%           t               vector specifying the time instances in seconds
%                           where the respective Fourier spectra were
%                           computed.
%           sideinfo.
%               stft.anaHop
%               stft.win
%               stft.zeroPad
%               stft.featureRate
%               stft.originalLength
%               stft.fsAudio
%               stft.fftShift
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
% check parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<3
    sideinfo=[];
end
if nargin<2
    parameter=[];
end
if nargin<1
    error('Please specify input data x.');
end

if ~isfield(parameter,'anaHop')
    parameter.anaHop = 2048;
end
if ~isfield(parameter,'win')
    parameter.win = win(4096,2); % hann window
end
if ~isfield(parameter,'zeroPad')
    parameter.zeroPad = 0;
end
if ~isfield(parameter,'fsAudio')
    parameter.fsAudio = 22050;
end
if ~isfield(parameter,'numOfFrames')
    parameter.numOfFrames = -1;
end
if ~isfield(parameter,'fftShift')
    parameter.fftShift = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% some pre calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% zeropad the window
w = parameter.win;
w = w(:);
zp = parameter.zeroPad;
w = [zeros(floor(zp/2),1);w;zeros(floor(zp/2),1)];
winLen = length(w);
winLenHalf = round(winLen/2);

fsAudio = parameter.fsAudio;
signalLength = length(x);
anaHop = parameter.anaHop;

% Pad the audio to center the windows and to avoid problems at the end
maxAnaHop = max(anaHop);
xPadded = [zeros(winLenHalf,1);x;zeros(winLen+maxAnaHop,1)];

% in case anaHop is a scalar, sample the window positions evenly in the
% input signal
if isscalar(anaHop)
    if parameter.numOfFrames >= 0
        numOfFrames = parameter.numOfFrames;
    else
        numOfFrames = floor((length(xPadded) - winLen)/anaHop + 1);
    end
    winPos = (0:numOfFrames-1) * anaHop + 1;
else
    if parameter.numOfFrames >= 0
        numOfFrames = parameter.numOfFrames;
    else
        numOfFrames = length(anaHop);
    end
    winPos = anaHop(1:numOfFrames);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spectrogram calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spec = zeros(winLenHalf+1,numOfFrames);
for i = 1 : numOfFrames
    xi = xPadded(winPos(i):winPos(i) + winLen - 1) .* w;
    if parameter.fftShift == 1
        xi = fftshift(xi);
    end
    Xi = fft(xi);
    spec(:,i) = Xi(1:winLenHalf+1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% axis calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = (winPos - 1) ./ fsAudio;
f = (0 : winLenHalf) * fsAudio / winLen;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update sideinfo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sideinfo.stft.anaHop = parameter.anaHop;
sideinfo.stft.win = parameter.win;
sideinfo.stft.zeroPad = parameter.zeroPad;
if isscalar(anaHop)
    sideinfo.stft.featureRate = fsAudio/anaHop;
else
    sideinfo.stft.featureRate = 'variable';
end
sideinfo.stft.originalLength = signalLength;
sideinfo.stft.fsAudio = parameter.fsAudio;
sideinfo.stft.fftShift = parameter.fftShift;

end
