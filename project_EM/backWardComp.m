function [lfpDat] = backWardComp(lfpDat)

if isfield(lfpDat,'dsFs')
    lfpDat.Fs = lfpDat.dsFs;
    lfpDat.trlTime = lfpDat.dsTrlTime;
    if min(lfpDat.trlTime) == -7
        for it = 1:length(lfpDat.LFPseg)
            lfpDat.LFPseg{it} = lfpDat.LFPseg{it}(lfpDat.trlTime>=-5,1:length(lfpDat.trlENC));
        end;
        lfpDat.trlTime = lfpDat.trlTime(lfpDat.trlTime>=-5);
    end;
end;