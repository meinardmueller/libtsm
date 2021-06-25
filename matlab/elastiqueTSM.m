function [y,sideinfo] = elastiqueTSM(x,s,parameter,sideinfo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: elastiqueTSM
% Date: 12-2013
% Programmer: Maximilian Schaefer, Jonathan Driedger
% http://www.audiolabs-erlangen.de/resources/MIR/TSMtoolbox/
%
% elastiqueTSM is a wrapper function which calls a online version of the
% proprietary TSM algorithm "elastique". For further information about the
% algorithm visit http://www.zplane.de/index.php?page=description-elastique
% and for further information about the API see the elastique_doku-folder
% at http://www.sonicapi.com.
%
% In order to be able to use the online service, you need a valid
% Access-Id which you can get by signing up at http://www.sonicapi.com.
% Afterwards you can insert your id in the code below.
%
% In addition you need the command-line tool cURL to transmit data in the
% network and here to upload your file to the server. cURL should be
% installed on every UNIX-based system. For Windows users, the curl.exe
% must be located inside the TSM Toolbox folder (we use it in version
% 7.33.0, 32bit).
% To Download cURL visit http://curl.haxx.se/download.html.
%
% Input:  x                 input signal.
%         s                 constant scaling factor. Due to limitations of
%                           the online service it needs to be in the range
%                           [0.25 .. 4]
%
%         parameter.
%             pitchShift    additionally also pitch shifts the signal by
%                           [-24 .. 24] semitones.
%             formantShift  additionally also pitch shifts the formants of
%                           by [-12 .. 12] semitones.
%             fsAudio       the sampling rate of the input audio signal x.
%
%
% Output: y                 the time-scale modified output signal.
%         sideinfo
%           elastique.pitchShift
%           elastique.formantShift
%           elastique.fsAudio
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
% access id
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
access_id = ''; % insert your id here

if isempty(access_id)
    warning('elastiqueTSM:noAccessId',...
        'You need to specify your access id.\nReturn y = 0.');
    y = 0;
    return;
end

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

if ~isfield(parameter,'pitchShift')
   parameter.pitchShift = 0;
end
if ~isfield(parameter,'formantShift')
   parameter.formantShift = 0;
end
if ~isfield(parameter,'fsAudio')
   parameter.fsAudio = 22050;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% some pre calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = 1/s; % invert s due to the input standard of the web service
pitchShift = parameter.pitchShift;
formantShift = parameter.formantShift;

% check parameter range
if s > 4 || s < 0.25
    s = 1;
    warning('elastiqueTSM:parameterOutOfRange',...
        'Stretching factor out of range, set to default: 1');
end
if parameter.pitchShift < -24 || parameter.pitchShift > 24
   parameter.pitchShift = 0;
   warning('elastiqueTSM:parameterOutOfRange',...
       'pitchShift out of range, set to default: 0');
end
if parameter.formantShift < -12 || parameter.formantShift > 12
   parameter.formantShift = 0;
   warning('elastiqueTSM:parameterOutOfRange',...
       'formant out of range, set to default: 0');
end

% write the file temporary to the workspace for upload
audiowrite('input.wav',x,parameter.fsAudio);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Upload the file to the elastique server and get the FILE_ID for process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% localize cURL on Windows machines (curl.exe must be placed inside of the
% TSM Toolbox folder)
if ispc
    fullFilenameMfile = mfilename('fullpath');
    locationToolbox = fileparts(fullFilenameMfile);
    locationToolbox = [locationToolbox '/'];
else
    locationToolbox = '';
end

% upload the file to server
cmd = [locationToolbox ...
    'curl -k "https://api.sonicapi.com/file/upload?access_id=',...
    access_id,'" -Ffile=@','input.wav'];
[~,w] = system(cmd);

% get the File-ID for the download
st = strfind(w,'file_id="') + length('file_id="');
en = strfind(w,'" status') - 1;


% if unsuccessfull
if isempty(st:en)
    warning('elastiqueTSM:curlError',...
        'It seems that cURL produced some error.\n%s\nReturn y = 0.',...
        w);
    y = 0;
    return;
end

file_id = w(st:en);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the process and download_link
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
serv = 'https://api.sonicapi.com/process/elastique?access_id=';
in = '&input_file=';
args = ['&tempo_factor=',num2str(s),'&pitch_semitones=',...
    num2str(pitchShift),'&formant_semitones=',...
    num2str(formantShift),'&format=wav'];
cmd = [serv,access_id,in,file_id,args];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process and download the file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[pathTmp,status] = urlwrite(cmd,'tmp.wav');
if status == 0
    warning('elastiqueTSM:elastiqueWebService',...
        'urlwrite produced some error.\nReturn y = 0.');
    y = 0;
    removeFiles();
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create the return value and delete tmp-file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% resampling
[y, fs] = audioread(pathTmp);
y = resample (y,parameter.fsAudio,fs,100);

removeFiles();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update sideinfo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sideinfo.elastique.pitchShift = parameter.pitchShift;
sideinfo.elastique.formantShift = parameter.formantShift;
sideinfo.elastique.fsAudio = parameter.fsAudio;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% remove temporary files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function removeFiles()
if(ispc == 1)
    system('del /F tmp.wav');
    system('del /F input.wav');
else
    system('rm -rf tmp.wav');
    system('rm -rf input.wav');
end
end
