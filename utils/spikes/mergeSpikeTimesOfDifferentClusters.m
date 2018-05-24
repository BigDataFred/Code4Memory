function [spkResp,bfLab] = mergeSpikeTimesOfDifferentClusters(sortedSpikes,trlSel,chanLab,selInfo,selIx)

spkResp = {};
bfLab = {};
it = 0;
cnt = 0;
while it <=length( selIx ) 
    
    it = it+1;
    if it > length( selIx )
        break;
    end;
    
    ts = [];    trl = [];	[wvf] = [];     cID = [];
    while true
        [spkDat] = sortedSpikes{selInfo(selIx(it),1)};
        [ts] = [ts spkDat.SpikeTimesSeg(spkDat.assignedClusterSeg == selInfo(selIx(it),2))];
        [trl] = [trl spkDat.trl(spkDat.assignedClusterSeg == selInfo(selIx(it),2))];
        [wvf] = [wvf;spkDat.wavf(spkDat.oriIx,:)];
        [cID] = [cID unique(spkDat.assignedClusterSeg(spkDat.assignedClusterSeg == selInfo(selIx(it),2)==1))];
        if it+1 > length( selIx )
            break;
        end;
        if  (selInfo(selIx(it),3) ~= 0) && (selInfo(selIx(it+1),3) ~= 0)
            it = it+1;
        else
            break;
        end;
    end;
    
    cnt = cnt+1;
    spkResp{cnt}.wvf = wvf(ismember(spkDat.assignedClusterSeg,cID),:);    
    for jt = 1:length(trlSel)
        spkResp{cnt}.spikeTimes{jt} = ts(trl==trlSel(jt)).*1e3;
    end;
    bfLab(cnt) = chanLab(selInfo(selIx(it),1));
    
end;