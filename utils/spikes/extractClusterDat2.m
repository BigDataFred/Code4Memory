function [cluDat] = extractClusterDat2(spkDat)
cnt = 0;
ts = {};
trl = {};
chanLab = {};
for it = 1:length(spkDat.sortedSpikes)
    dat = spkDat.sortedSpikes{it};
    
    uID  = unique(dat.assignedClusterSeg);
    
    for jt = 1:length( uID)
        cnt = cnt+1;
        trl{cnt} = dat.trl(dat.assignedClusterSeg == uID(jt));
        ts{cnt} = dat.SpikeTimesSeg(dat.assignedClusterSeg == uID(jt)).*1e3;                
        chanLab{cnt} = spkDat.chanLab{it};
    end;
end;
cluDat.trl = trl;
cluDat.ts = ts;
cluDat.chanLab = chanLab;
return;