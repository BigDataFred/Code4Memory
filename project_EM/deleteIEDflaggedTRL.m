function [trlENC,delIx,MicTrlIx1,MicTrlIx2] = deleteIEDflaggedTRL(micDat,ixIEDdat,macDat,trlENC)

dum = micDat{1}.save_data{1}{1}{1}.trial;
ntrl = length(dum);
% if length(trlENC)*2~=ntrl
%     error('inconsistent trial assignment detected');
% end;

%% delete those trial labels that have been flaged by IED detection
% search for indexes flaged with IED
delIx = [];
if isfield(micDat{1},'delIx')
    for it = 1:length(micDat)% loop over the MWs
        delIx = [delIx micDat{it}.delIx];
    end;
elseif ~isempty(ixIEDdat) % else use average MW-LFP
    delIx = ixIEDdat.delIx;
end;

% get bad trials from macro-LFP-data
if ~isempty( macDat )
    delIx = [delIx macDat.delIx'];
end;

%delete Encoding trial labels that have been flaged up
delIx = unique(delIx);
trlENC(ismember(trlENC,delIx))=[];

% extract index of Encoding trials that correspond to clean MW data
selIx = [];
if isfield(micDat{1},'selIx')
    for it = 1:length(micDat)% loop over the MWs
        selIx = [selIx micDat{it}.selIx];
    end;
    selIx = unique(selIx);
elseif isstruct(ixIEDdat) && ~isempty(ixIEDdat.selIx)% use avg-MW
    selIx = ixIEDdat.selIx;
else
    selIx = [];
end;

if isempty(selIx)
    [MicTrlIx1] = trlENC;
    [MicTrlIx2] = find(ismember(1:length(micDat{1}.save_data{1}{1}{1}.trial),MicTrlIx1));
else
    [MicTrlIx1] = intersect(trlENC,selIx);
    [MicTrlIx2] = find(ismember(selIx,MicTrlIx1));
end;