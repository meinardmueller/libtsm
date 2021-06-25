function [y,sideinfo] = pvTSM(x,s,parameter,sideinfo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: pvTSM
% Date: 03-2014
% Programmer: Jonathan Driedger
% http://www.audiolabs-erlangen.de/resources/MIR/TSMtoolbox/
%
% The phase vocoder is a time-scale modification algorithm. It rescales the
% time-axis of the input signal x according to the time-stretch function s
% without altering the pitch of x. The phase vocoder is particulary well
% suited for input signals of harmonic nature. For further details see the
% paper "Phase Vocoder" by Flanagan and Golden. Improvements to the basic
% phase vocoder technique, in particular a technique called 'identity phase
% locking' has been proposed in "Improved Phase Vocoder Time-Scale
% Modification of Audio" by Laroche and Dolson and can optionally also be
% used in this implementation.
%
% Input:  x                 input signal.
%         s                 time-stretch function. Either a constant
%                           scaling factor or a n x 2 matrix representing a
%                           set of n anchorpoints relating sample positions
%                           in the input signal with sample positions in
%                           the output signal.
%         parameter.
%            synHop         hop size of the synthesis window.
%            win            the analysis and synthesis window for the stft.
%            zeroPad        number of zeros that should be padded to the
%                           window to increase the fft size and therefore
%                           the frequency resolution.
%            restoreEnergy  set to 1 in case the istft should account for
%                           a potential energy loss of the output signal.
%            fftShift       set to 1 in case the stft and istft should
%                           apply a circular shift by half the frame length
%                           to each frame prior to their application.
%            phaseLocking   set to 1 in case identity phase locking should
%                           be applied.
%
% Output: y                 the time-scale modified output signal.
%         sideinfo.
%             pvTSM.synHop
%             pvTSM.win
%             pvTSM.zeroPad
%             pvTSM.restoreEnergy
%             pvTSM.fftShift
%             pvTSM.phaseLocking
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
if nargin < 4
    sideinfo = [];
end
if nargin < 3
    parameter = [];
end
if nargin<2
    error('Please specify input data x and s.');
end

if ~isfield(parameter,'synHop')
    parameter.synHop = 512;
end
if ~isfield(parameter,'win')
    parameter.win = win(2048,1); % sin window
end
if ~isfield(parameter,'zeroPad')
    parameter.zeroPad = 0;
end
if ~isfield(parameter,'restoreEnergy')
    parameter.restoreEnergy = 0;
end
if ~isfield(parameter,'fftShift')
    parameter.fftShift = 0;
end
if ~isfield(parameter,'phaseLocking')
    parameter.phaseLocking = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% some pre calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
synHop = parameter.synHop;
w = parameter.win;
zp = parameter.zeroPad;
w = [zeros(floor(zp/2),1);w;zeros(floor(zp/2),1)];
winLen = length(w);
winLenHalf = round(winLen/2);
phaseLocking = parameter.phaseLocking;
numOfChan = size(x,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% time-stretch function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isscalar(s)
    anchorPoints = [1 1; size(x,1) ceil(s*size(x,1))];
else
    anchorPoints = s;
end
outputLength = anchorPoints(end,2);
synWinPos = 1:synHop:outputLength + ... % Positions of the synthesis
    winLenHalf;                         % windows in the output
anaWinPos = ...                         % Positions of the analysis windows
    round(interp1(anchorPoints(:,2),... % in the input
    anchorPoints(:,1),synWinPos,...
    'linear','extrap'));
anaHop = [0 anaWinPos(2:end)-...        % Analysis hopsizes
    anaWinPos(1:end-1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% phase vocoder
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y = zeros(outputLength,numOfChan);      % Initialize output

for c = 1 : numOfChan                   % loop over channels
xC = x(:,c);

% stft
parStft.anaHop = anaWinPos;
parStft.win = parameter.win;
parStft.zeroPad = parameter.zeroPad;
parStft.fftShift = parameter.fftShift;
X = stft(xC,parStft);

% phase adaption
Y = zeros(size(X));                     % The spectrogram of the output
Y(:,1) = X(:,1);                        % Phase initialization

N = length(w);
k = (0:N/2)';                           % Center frequencies of the N/2+1
                                        % first bins of the spectrum
                                        % in 'oscillations per frame'
omega = 2*pi*k/N;                       % Phase advances per sample for the
                                        % frequencies k
for i = 2 : size(X,2);
    dphi = omega * anaHop(i);           % Expected phase advances from the
                                        % last to the current input frame
    phCurr = angle(X(:,i  ));           % Phases of the current input frame
    phLast = angle(X(:,i-1));           % Phases of the last input frame

    hpi = (phCurr - phLast) - dphi;     % Heterodyned phase increments
    hpi = hpi - 2 * pi * ...
        round(hpi/(2*pi));              % Reduce to the range -pi:pi

    ipa_sample = (omega+hpi/anaHop(i)); % Instantaneous phase advances per
                                        % sample
    ipa_hop = ipa_sample * synHop;      % Instantaneous phase advances per
                                        % synthesis hopsize

    phSyn = angle(Y(:,i-1));            % Phases of the last synthesized
                                        % frame

    % We now compute a phasor that rotates the phase angles of the current
    % input frame by angles theta such that no phase discontinuities occur
    % when resynthesizing the resulting spectrogram with the synthesis
    % hopsize
    if phaseLocking == 0                % The standard phase vocoder: the
                                        % phase continuity of every bin is
                                        % preserved separately

        theta = phSyn + ...             % Phases of the last output frame
            ipa_hop - ...               % Instantaneous phase advance
            phCurr;                     % Phases of the current input frame
        phasor = exp(1i*theta);

    else                                % Phase vocoder with identity phase
                                        % locking: the phase relationships
                                        % from the input frame are
                                        % partially preserved by 'locking'
                                        % the phases of bins in the region
                                        % of influence of a peak in the
                                        % sprectrum to the phase of the
                                        % peak's bin

        [p,irS,irE] = findPeaks(X(:,i));% Get the peaks in the spectrum
                                        % together with their regions of
                                        % influence

        theta = zeros(size(Y(:,i)));
        for n = 1 : length(p)
            theta(irS(n):irE(n)) = ...
               phSyn(p(n)) + ...        % Phases of the last output frame
               ipa_hop(p(n)) - ...      % Instantaneous phase advance
               phCurr(p(n));            % Phases of the current input frame
        end
        phasor = exp(1i*theta);

    end
    Y(:,i) = phasor .* X(:,i);          % Rotate angles of the current
                                        % input frame to get the output
                                        % frame
end

% istft
parIstft.synHop = synHop;
parIstft.win = parameter.win;
parIstft.zeroPad = parameter.zeroPad;
parIstft.numOfIter = 1;
parIstft.origSigLen = outputLength;
parIstft.restoreEnergy = parameter.restoreEnergy;
parIstft.fftShift = parameter.fftShift;
yC = istft(Y,parIstft);

y(:,c) = yC;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update sideinfo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sideinfo.pvTSM.synHop = parameter.synHop;
sideinfo.pvTSM.win = parameter.win;
sideinfo.pvTSM.zeroPad = parameter.zeroPad;
sideinfo.pvTSM.restoreEnergy = parameter.restoreEnergy;
sideinfo.pvTSM.fftShift = parameter.fftShift;
sideinfo.pvTSM.phaseLocking = parameter.phaseLocking;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find peaks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [peaks,inflRegionStart,inflRegionEnd] = findPeaks(spec)
% An index in spec is considered a peak if its value is the largest among
% its four nearest neighbours

magSpec = abs(spec);
magSpec = magSpec(:);
magSpecPadded = [zeros(2,1);magSpec;zeros(2,1)];

peaks = find(magSpecPadded(5:end  )  <  magSpecPadded(3:end-2) &...
             magSpecPadded(4:end-1)  <  magSpecPadded(3:end-2) &...
             magSpecPadded(2:end-3)  <  magSpecPadded(3:end-2) &...
             magSpecPadded(1:end-4)  <  magSpecPadded(3:end-2));

inflRegionStart = zeros(size(peaks));
inflRegionEnd   = zeros(size(peaks));

if isempty(peaks)
    return;
end

inflRegionStart(1) = 1;
inflRegionStart(2:end) = ceil((peaks(2:end) + peaks(1:end-1))/2);
inflRegionEnd(1:end-1) = inflRegionStart(2:end) - 1;
inflRegionEnd(end) = length(inflRegionEnd);

end
