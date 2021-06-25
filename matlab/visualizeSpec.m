function visualizeSpec(spec,parameter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: visualizeSpec
% Date: 03-2014
% Programmer: Jonathan Driedger
% http://www.audiolabs-erlangen.de/resources/MIR/TSMtoolbox/
%
% Visualizes a waveform.
%
% Input:    spec             spectrogram.
%           parameter.
%               fAxis        frequency axis of the spectrogram.
%               tAxis        timeAxis of the spectrogram.
%               logComp      factor for logarithmic compression of the data
%                            for a clearer visualization.
%               visFreqRange frequency range that should be visualized.
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
if nargin<2
    parameter=[];
end
if nargin<1
    error('Please specify input data x.');
end

if ~isfield(parameter,'fAxis')
    parameter.fAxis = 1:size(spec,1);
end
if ~isfield(parameter,'tAxis')
    parameter.tAxis = 1:size(spec,2);
end
if ~isfield(parameter,'logComp')
    parameter.logComp = 10;
end
if ~isfield(parameter,'visFreqRange')
    parameter.visFreqRange = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = parameter.fAxis;
t = parameter.tAxis;
logComp = parameter.logComp;
visFreqRange = parameter.visFreqRange;

figure;
g = gray;
g = g(end:-1:1,:);

imagesc(t,f,log(1 + logComp*abs(spec)));
colormap(g);
axis xy;

if ~isempty(visFreqRange)
    ylim(visFreqRange)
end

xlabel('Time [sec]');
ylabel('Frequency [Hz]');

end
