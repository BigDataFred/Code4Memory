function [LFPr] = computeBFcorrLFP_datModeA(MWdat)

%% extract the median LFP
Fs = 1/1e-3;
LFPr = [];
for kt = 1:length(MWdat.trial)
    nC = length(MWdat.label);
    for it = 1:length(MWdat.label)
        x = MWdat.trial{kt}(it,:)';
        idx = 1:nC;
        idx = setdiff(idx,it);
        for jt = 1:nC%length(idx)
            y = MWdat.trial{kt}(jt,:)';%idx(jt)
            LFPr(kt,it,jt) = corr(x,y);
        end;
    end;
end;