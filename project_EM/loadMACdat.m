function [MACdat] = loadMACdat(micDat,macDat,trlENC)
% list of trials that were retained for macro-data
macIx = setdiff(1:length([1:length(macDat.macroLFPdat.trial) macDat.delIx']),macDat.delIx);

% list of trials in macro data
macIx2 = 1:length( macDat.macroLFPdat.trial );

% extrat trials from macro-data corresponding to Encoding
dum = intersect(trlENC,macIx);
[selIx] = find(ismember(macIx,dum));

% extract trials from macro-data that correspond to trials rejected in
% micro data
delIx = [];
for it = 1:length(micDat)
    delIx = [delIx micDat{it}.delIx];
end;
delIx = unique(delIx);

dum = intersect(macIx,delIx);
[delIx] = find(ismember(macIx,dum));

% readout macro data trials corresponding to Encoding
MacTrlIx = macIx2(setdiff(selIx,delIx));

%%
cfg                     = [];
cfg.trials              = MacTrlIx;

[MACdat] = ft_selectdata( cfg , macDat.macroLFPdat );