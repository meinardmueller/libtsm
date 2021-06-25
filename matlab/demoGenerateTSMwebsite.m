%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: demoGenerateTSMWebsite
% Date: 12-2013
% Programmer: Jonathan Driedger
% http://www.audiolabs-erlangen.de/resources/MIR/TSMtoolbox/
%
% This is a demo script which shows how the 'TSM toolbox' can be used to
% generate a website for comparing the TSM results of the various TSM
% algorithms included in the toolbox.
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
pathData = '../data/';                        % path to the source files
pathOutput = '../output/';
pathHtml = [pathOutput 'website/'];        % path to the website folder
nameTable = 'tables.html';
nameIndexTemplate = 'indexTemplate.html';
nameStyle = 'style.css';
nameFiles = {                              % files that should be processed
    'Bongo.wav';
    'CastanetsViolin.wav';
    'DrumSolo.wav';
    'Glockenspiel.wav';
    'Stepdad.wav';
    'Jazz.wav';
    'Pop.wav';
    'SingingVoice.wav';
    'SynthMono.wav';
    'SynthPoly.wav';
    };

stretchingFactors = [0.5 1.2 1.8];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% file handling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate website directory
if ~isdir(pathHtml)
    mkdir(pathHtml)
end

% html index template
fid = fopen([pathData nameIndexTemplate],'r');
index_template = fread(fid, '*char')';
fclose( fid );

% html index file where the code for the table will be put in
index_out = fopen([ pathHtml 'index.html'], 'w');

% create the html file for the tables
hTables = fopen([pathHtml nameTable],'w+t');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate the tables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
warning('OFF','MATLAB:audiovideo:audiowrite:dataClipped');

for s = stretchingFactors
    disp(['Generating table and audio files for stretching factor ' ...
        num2str(s)]);

    % write title of the table
    fprintf(hTables,'<h3>Constant stretching factor of &alpha;=%s</h3>',...
        num2str(s));

    % write the head of the table
    fprintf(hTables,'<center>\n');
    fprintf(hTables, '<table>\n');

    % write the head row of the table
    fprintf(hTables, '<tr>\n');
    fprintf(hTables, '<th><b>Name</b></th>\n');
    fprintf(hTables, '<th><b>Original</b></th>\n');
    fprintf(hTables, '<th><b>OLA</b></th>\n');
    fprintf(hTables, '<th><b>WSOLA</b></th>\n');
    fprintf(hTables, '<th><b>Phase Vocoder</b></th>\n');
    fprintf(hTables, '<th><b>Phase Vocoder <br>(phase locking)</b></th>\n');
    fprintf(hTables, '<th><b>TSM based on <br>HPSS</b></th>\n');
    fprintf(hTables, '<th><b>&eacute;lastique</b></th>\n');
    fprintf(hTables, '</tr>\n');

    % separator line
    fprintf(hTables, '<tr align="center">\n');
    fprintf(hTables, '<th class="separator"></th>\n');
    fprintf(hTables, '<th class="separator"></th>\n');
    fprintf(hTables, '<th class="separator"></th>\n');
    fprintf(hTables, '<th class="separator"></th>\n');
    fprintf(hTables, '<th class="separator"></th>\n');
    fprintf(hTables, '<th class="separator"></th>\n');
    fprintf(hTables, '<th class="separator"></th>\n');
    fprintf(hTables, '<th class="separator"></th>\n');
    fprintf(hTables, '</tr>\n');

    % iterate over all items
    for i = 1: length(nameFiles)
        name = nameFiles{i};
        disp(['  ' name]);

        % construct filenames
        name_ORIG = [name(1:end-4) '_ORIG.wav'];
        name_OLA = [name(1:end-4) '_' sprintf('%0.2f',s) '_OLA.wav'];
        name_WSOLA = [name(1:end-4) '_' sprintf('%0.2f',s) '_WSOLA.wav'];
        name_PV = [name(1:end-4) '_' sprintf('%0.2f',s) '_PV.wav'];
        name_PVpl = [name(1:end-4) '_' sprintf('%0.2f',s) '_PVpl.wav'];
        name_HP = [name(1:end-4) '_' sprintf('%0.2f',s) '_HP.wav'];
        name_ELAST = [name(1:end-4) '_' sprintf('%0.2f',s) '_ELAST.wav'];

        % generate audio files
        % Original
        [x,fsAudio] = audioread([pathData name]);
        audiowrite([pathHtml name_ORIG],x,fsAudio);

        % OLA
        paramOLA.tolerance = 0;
        paramOLA.synHop = 128;
        paramOLA.win = win(256,2); % hann window
        yOLA = wsolaTSM(x,s,paramOLA);
        audiowrite([pathHtml name_OLA],yOLA,fsAudio);

        % WSOLA
        yWSOLA = wsolaTSM(x,s);
        audiowrite([pathHtml name_WSOLA],yWSOLA,fsAudio);

        % Phase Vocoder
        yPV = pvTSM(x,s);
        audiowrite([pathHtml name_PV],yPV,fsAudio);

        % Phase Vocoder with identity phase locking
        paramPVpl.phaseLocking = 1;
        yPVpl = pvTSM(x,s,paramPVpl);
        audiowrite([pathHtml name_PVpl],yPVpl,fsAudio);

        % TSM based on HPSS
        yHP = hpTSM(x,s);
        audiowrite([pathHtml name_HP],yHP,fsAudio);

        % elastique
        yELAST = elastiqueTSM(x,s);
        audiowrite([pathHtml name_ELAST],yELAST,fsAudio);


        % generate table entries
        % setup new table row
        fprintf(hTables, '<tr align="center">\n');
        fprintf(hTables, '<td align="left" class="name">%s</td>\n',...
            name(1:end-4));

        % Original
        fprintf(hTables,'<td><a href="%s">[wav]</a></td>\n',name_ORIG);

        % OLA
        fprintf(hTables,'<td><a href="%s">[wav]</a></td>\n',name_OLA);

        % WSOLA
        fprintf(hTables,'<td><a href="%s">[wav]</a></td>\n',name_WSOLA);

        % Phase Vocoder
        fprintf(hTables,'<td><a href="%s">[wav]</a></td>\n',name_PV);

        % Phase Vocoder with identity phase locking
        fprintf(hTables,'<td><a href="%s">[wav]</a></td>\n',name_PVpl);

        % TSM based on HPSS
        fprintf(hTables,'<td><a href="%s">[wav]</a></td>\n',name_HP);

        % elastique
        fprintf(hTables,'<td><a href="%s">[wav]</a></td>\n',name_ELAST);

        % close table row
        fprintf(hTables, '</tr>\n');

    end

    % closing the table
    fprintf(hTables, '</table>\n');
    fprintf(hTables, '</center>\n');

end

fclose(hTables);


% load the written table
fid = fopen([pathHtml nameTable],'r');
table = fread(fid, '*char')';
fclose( fid );

% replace the labels in indexTemplate
index_template = regexprep( index_template, 'LABEL_TABLE', table);

% write the index.html
fprintf(index_out,'%s',index_template);
fclose(index_out);

% copy the stylesheet to the html folder
copyfile([pathData nameStyle],[pathHtml nameStyle]);
