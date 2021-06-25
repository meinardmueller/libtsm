function [y,sideinfo] = pvIntTSM(x,s,parameter,sideinfo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: pvIntTSM
% Date: 05-2015
% Programmer: Jonathan Driedger
%
% The phase vocoder is a time-scale modification algorithm. It rescales the
% time-axis of the input signal x according to the time-stretch function s
% without altering the pitch of x. The phase vocoder is particulary well
% suited for input signals of harmonic nature. For further details see
% [1966_FlanaganGolden_PhaseVocoder_BellTechnicalJournal].
%
% Input:    x               single-channel signal.
%           s               time-stretch function. Either a constant
%                           scaling factor or a n x 2 matrix representing a
%                           set of n anchorpoints relating sample positions
%                           in the input signal with sample positions in
%                           the output signal.
%           parameter.
%               synHop      hop size of the synthesis window.
%               window      the analysis and synthesis window for the stft.
%               zeroPad     number of zeros that should be padded to the
%                           window to increase the fft size and therefore
%                           the frequency resolution.
%
% Output:   y               the time-scale modified output signal.
%           sideinfo.
%               pvTsm.synHop
%               pvTsm.window
%               pvTsm.zeroPad
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

if ~isfield(parameter,'synHop')
    parameter.synHop = 512;
end
if ~isfield(parameter,'win')
    parameter.win = win(2048,2); % hann window
end
if ~isfield(parameter,'zeroPad')
%     parameter.zeroPad = (s-1)*length(parameter.win);
    parameter.zeroPad = s*length(parameter.win)/2;
%      parameter.zeroPad = 0;
end
if ~isfield(parameter,'restoreEnergy')
    parameter.restoreEnergy = 0;
end
if ~isfield(parameter,'fftShift')
    parameter.fftShift = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% some pre calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
synHop = parameter.synHop;
w = parameter.win;
zp = parameter.zeroPad;
% zp = (s-1) * length(w);
w = [zeros(floor(zp/2),1);w;zeros(floor(zp/2),1)];
winLen = length(w);
winLenHalf = round(winLen/2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% time-stretch function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isscalar(s) && mod(s,1) == 0 && s>=1
    anchorPoints = [1 1; length(x) ceil(s*length(x))];
else
    error('pvIntTSM works just for integer stretching factors s >= 1.')
end

while mod(synHop,s) ~= 0
    synHop = synHop + 1;
end

outputLength = anchorPoints(end,2);
outputWindowPos = 1:synHop:outputLength + winLenHalf;
inputWindowPos = (outputWindowPos-1)./s + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parStft.anaHop = inputWindowPos;
parStft.win = parameter.win;
parStft.zeroPad = parameter.zeroPad;
parStft.fftShift = parameter.fftShift;
X = stft(x,parStft);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% phase adaption
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = abs(X) .* exp(1i * s * angle(X));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% istft
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parIstft.synHop = synHop;
parIstft.win = parameter.win;
parIstft.zeroPad = parameter.zeroPad;
parIstft.numberOfIter = 1;
parIstft.origSigLen = outputLength;
parIstft.restoreEnergy = parameter.restoreEnergy;
parIstft.fftShift = parameter.fftShift;
y = istft(Y,parIstft);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update sideinfo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sideinfo.pvTsm.synHop =  parameter.synHop;
sideinfo.pvTsm.win = parameter.win;
sideinfo.pvTsm.zeroPad = parameter.zeroPad;

end
