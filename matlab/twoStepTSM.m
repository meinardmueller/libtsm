function [y,sideinfo] = twoStepTSM(x,s,parameter,sideinfo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: twoStepTSM
% Date: 05-2015
% Programmer: Jonathan Driedger
%
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

if ~isfield(parameter,'roughTSM')
    parameter.roughTSM = @pvIntTSM;
end
if ~isfield(parameter,'exactTSM')
    parameter.exactTSM = @hpTSM;
end
if ~isfield(parameter,'cascOrder')
    parameter.cascOrder = 'exactCoarse';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TSM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sRough = max(1,round(s));
sExact = s/sRough;

roughTSM = parameter.roughTSM;
exactTSM = parameter.exactTSM;

switch parameter.cascOrder

    case 'exactCoarse'
        yTmp = exactTSM(x,sExact);
        y = roughTSM(yTmp,sRough);

    case 'coarseExact'
        yTmp = roughTSM(x,sRough);
        y = exactTSM(yTmp,sExact);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update sideinfo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sideinfo.twoStepTsm.yTmp =  yTmp;
sideinfo.twoStepTsm.sRough = sRough;
sideinfo.twoStepTsm.sExact = sExact;

end
