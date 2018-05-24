function [macro_filename] = find_macro_filename_EM(mic_ts,pID)
%
% This function converts the NLX timestamps into EDF format and performs a
% lookup in the list of EDF files for a corresponding match. A filename is
% returned in the output that contains macro-channel data, which corresponds
% with the timerange of the NLX timestamps.
%
% Input:    1. mic_ts: Neuralynx timestamps, string format must be yyyy-mm-dd_hh-mm-ss
%           2. pId: dataset ID, string format must be 'Pxx' with xx ranging from 01-n
%
% Output:   1. macro_filename: the filename of the .edf file that matches with the Nlx timestamps


%F. Roux aka FRO, 
%University of Birmingham 2016

%% convert NLX timestamps to correspond with EDF header format
ix = regexp(mic_ts,'_');

d = mic_ts(1:ix-1);% extract the date (yyyy-mm-dd)
h = mic_ts(ix+1:length(mic_ts));% extract the time (hh-mm-ss)

% parse the string with the date info
ix = regexp(d,'[-]');
[pVal1] = parse_str2val(d,ix);

% parse the string with the time info
ix = regexp(h,'[-]');
[pVal2] = parse_str2val(h,ix);

% adapte the NLX timestamps to EDF header format
d = [pVal1{2},'/',pVal1{3},'/',pVal1{1}];
h = [pVal2{1},':',pVal2{2},':',pVal2{3}];

%% in the list of .edf files, search for an entry that matches NLX timestamps
[macro_filename] = table_lookup(d,h,pID);

%% 
function [parsedVal] = parse_str2val(str,ix)
% parse string into discrete elements
dum = [];
for it = 1:length(ix)+1    
    if (it >1) && (it <=length(ix))
        dum = [dum str(ix(it-1)+1:ix(it)-1)];
    elseif (it ==1)
        dum = [dum {str(1:ix(it)-1)}];
    else
        dum = [dum str(ix(it-1)+1:length(str))];
    end;    
end;
parsedVal = dum;
for it = 1:length(parsedVal)
    parsedVal(it) = {num2str(str2double(parsedVal{it}))};
end;

%%
function [filename] = table_lookup(d,h,pID)

p2lookup = ['/media/rouxf/rds-share/Archive/MACRO/',pID,'/EDF/'];% path where lists for lookup are located
lookupFile = [pID,'_macro_micro_lookup.txt'];% this files contains a list of all EDF files
missingFile = [pID,'_missing_macro_files.txt'];% list of NLX timestamps for which no corresponding EDF file was found

%check if file exist
lookupFile = dir([p2lookup,lookupFile]);
missingFile = dir([p2lookup,missingFile]);

% read the file info
fid_ = fopen([p2lookup,lookupFile.name],'r');
fidx = fopen([p2lookup,missingFile.name],'r');

% import the lookup table
lookupDat = textscan(fid_,'%s');
lookupDat = lookupDat{:};
lookupDat = reshape(lookupDat,[5 length(lookupDat)/5])';

for it = 1:size( lookupDat(:,2) ,1)
    ix = regexp(lookupDat{it,2},'[.]');
    lookupDat(it,2) = {lookupDat{it,2}(1:ix-1)};
end;

% import the missing-data list
missingDat = textscan(fidx,'%s');
missingDat = missingDat{:};
missingDat = reshape(missingDat,[3 length(missingDat)/3])';

for it = 1:size( missingDat(:,2) ,1)
    ix = regexp(missingDat{it,2},'[.]');
    missingDat(it,2) = {missingDat{it,2}(1:ix-1)};
end;

% search for entries in the loopup table that overlab with the date and
% time of the NLX timestamps
[ix] = intersect(find(strcmp(d,lookupDat(:,1))),find(strcmp(h,lookupDat(:,2))));
if ~isempty(ix)
    [filename] = lookupDat{ix,end};
else
    %[ix] = intersect(find(strcmp(d,missingDat(:,1))),find(strcmp(h,missingDat(:,2))));
    [filename] = [];
end;

fclose('all');
return;