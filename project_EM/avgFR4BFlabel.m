function [avgFR,semFR,BFid] = avgFR4BFlabel(spkDat,BFlab)

[BFid] = unique(BFlab);

[avgFR] = cell(length(BFid),1);
[semFR] = cell(length(BFid),1);

for curBF = 1:length(BFid)
    ix = find(strcmp(BFlab,BFid(curBF)));
    avgFR{curBF} = mean(spkDat.fr(ix,:),1);
    semFR{curBF} = std(spkDat.fr(ix,:),0,1)./sqrt(length(ix)-1);
end;

return;