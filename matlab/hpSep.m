function [xHarm,xPerc,sideinfo] = hpSep(x,parameter,sideinfo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: hpSep
% Date: 03-2014
% Programmer: Jonathan Driedger
% http://www.audiolabs-erlangen.de/resources/MIR/TSMtoolbox/
%
% Seperates a given audio signal into a harmonic and a percussive component
% according to the paper "Harmonic/Percussive Separation using Median
% Filtering" by Fitzgerald.
%
% Input:  x                 input signal.
%         parameter.
%          anaHop           the stft hop size of the analysis window.
%          win              the stft analysis window used for windowing the
%                           input signal.
%          zeroPad          number of zeros that should be padded to the
%                           window to increase the fft size and therefore
%                           the frequency resolution.
%          filLenHarm       length of the median filter in time direction.
%          filLenPerc       length of the median filter in frequency
%                           direction.
%          maskingMode      either 'binary' or 'relative'. Specifies if a
%                           binary or a relative weighting mask should be
%                           applied to the spectrogram to perform the
%                           separation.
%
% Output: xHarm             the harmonic component of the input signal x.
%         xPerc             the percussive component of the input signal x.
%
%         sideinfo.
%            hpSep.stftAnaHop
%            hpSep.win
%            hpSep.zeroPad
%            hpSep.filLenHarm
%            hpSep.filLenPerc
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
if nargin < 3
    sideinfo = [];
end
if nargin < 2
    parameter = [];
end
if nargin<2
    error('Please specify input data x.');
end

if ~isfield(parameter,'anaHop')
    parameter.anaHop = 256;
end
if ~isfield(parameter,'win')
    parameter.win = win(1024,2); % hann window
end
if ~isfield(parameter,'zeroPad')
    parameter.zeroPad = 0;
end
if ~isfield(parameter,'filLenHarm')
    parameter.filLenHarm = 10;
end
if ~isfield(parameter,'filLenPerc')
    parameter.filLenPerc = 10;
end
if ~isfield(parameter,'maskingMode')
    parameter.maskingMode = 'binary';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% some pre calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
anaHop = parameter.anaHop;
w = parameter.win;
zeroPad = parameter.zeroPad;
filLenHarm = parameter.filLenHarm;
filLenPerc = parameter.filLenPerc;
maskingMode = parameter.maskingMode;
numOfChan = size(x,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% harmonic-percussive separation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xHarm = zeros(size(x,1),numOfChan);      % Initialize output
xPerc = zeros(size(x,1),numOfChan);      % Initialize output

for c = 1 : numOfChan                   % loop over channels
xC = x(:,c);

% stft
parStft.anaHop = anaHop;
parStft.win = w;
parStft.zeroPad = zeroPad;
spec = stft(xC,parStft);
magSpec = abs(spec);

% harmonic-percussive separation
magSpecHarm = medianFilter(magSpec,filLenHarm,2);
magSpecPerc = medianFilter(magSpec,filLenPerc,1);

switch maskingMode
    case 'binary'
        maskHarm = magSpecHarm >  magSpecPerc;
        maskPerc = magSpecHarm <= magSpecPerc;

    case 'relative'
        maskHarm = magSpecHarm ./ (magSpecHarm + magSpecPerc + eps);
        maskPerc = magSpecPerc ./ (magSpecHarm + magSpecPerc + eps);

    otherwise
        error('maskingMode must either be "binary" or "relative"');
end

specHarm = maskHarm .* spec;
specPerc = maskPerc .* spec;

% istft
parIstft.synHop = parameter.anaHop;
parIstft.win = parameter.win;
parIstft.zeroPad = parameter.zeroPad;
parIstft.numOfIter = 1;
parIstft.origSigLen= length(x);
xHarmC = istft(specHarm,parIstft);
xPercC = istft(specPerc,parIstft);

xHarm(:,c) = xHarmC;
xPerc(:,c) = xPercC;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update sideinfo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sideinfo.hpSep.stftAnaHop = parameter.anaHop;
sideinfo.hpSep.win = parameter.win;
sideinfo.hpSep.zeroPad = parameter.zeroPad;
sideinfo.hpSep.filLenHarm = parameter.filLenHarm;
sideinfo.hpSep.filLenPerc = parameter.filLenPerc;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% median filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Y = medianFilter(X,len,dim)

s = size(X);
Y = zeros(s);

switch dim
    case 1
        XPadded = [zeros(floor(len/2),s(2));X;zeros(ceil(len/2),s(2))];
        for i = 1 : s(1)
            Y(i,:) = median(XPadded(i:i+len-1,:),1);
        end

    case 2
        XPadded = [zeros(s(1),floor(len/2)) X zeros(s(1),ceil(len/2))];
        for i = 1 : s(2)
            Y(:,i) = median(XPadded(:,i:i+len-1),2);
        end

    otherwise
        error('unvalid dim.')
end

end
