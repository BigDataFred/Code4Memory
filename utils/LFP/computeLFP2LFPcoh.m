function [LFP2LFPcoh] = computeLFP2LFPcoh(LFPdat,Fs)

TW = 5;

params                  = [];
params.Fs               = Fs;
params.pad              = 2;
params.tapers           = [TW 2*TW-1];
params.fpass            = [1 30];
params.trialave         = 1;
params.err              = 0;

C= cell(length(LFPdat),length(LFPdat)-1);
fx = [];
for it = 1:length(LFPdat)
    idx = setdiff(1:length(LFPdat),it);
    for jt = 1:length(idx)
        
        [C{it,jt},~,~,~,~,fx] = coherencyc(LFPdat{it},LFPdat{idx(jt)},params);
        
    end;
end;

LFP2LFPcoh.C = C;
LFP2LFPcoh.fx = fx;