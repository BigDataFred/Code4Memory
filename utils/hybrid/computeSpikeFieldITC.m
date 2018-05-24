function [ITC,freqAx,timeAx,pval] = computeSpikeFieldITC(data_spk,data_lfp, Fs,timeGrid,bootStrap)
%%

return;

return;

function [STA] = plugDat2FT(STA,timeGrid,Fs)

ntrl = size(STA,1);

dum = [];
dum.fsample = Fs;
dum.time = {};
dum.trial = {};
dum.label = {'STAchan'};

for it = 1:ntrl
    dum.time{it} = timeGrid;
    dum.trial{it} = STA(it,:);
end;
STA = dum;
