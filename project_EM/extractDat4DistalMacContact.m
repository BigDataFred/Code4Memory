function [macDat] = extractDat4DistalMacContact(macDat)

%% extract the labels of the most distal contacts on the Macro channels
chanNum = zeros(length(macDat.label),1);
chanID  = cell(length(macDat.label),1);

for it = 1:length(macDat.label)
    
    selIx = regexp(macDat.label{it},'\d');
    chanNum(it) = str2double(macDat.label{it}(selIx));
    chanID(it) = {macDat.label{it}(1:min(selIx)-1)};
end;

[macID] = unique(chanID);

selIx = zeros(length(macID),1);
for it = 1:length(macID)
    
    chanIx = find(strcmp(chanID,macID(it)));
    
    [~,ix] = min(chanNum(chanIx));
    
    selIx(it) = chanIx(ix);
    
end;

%%
cfg                     = [];
cfg.channel             = macDat.label(selIx);

[macDat] = ft_selectdata( cfg , macDat );
