function [avgFR,BFid] = avgFR4BFlabel(spkDat,BFlab)

[BFid] = unique(BFlab);

[avgFR] = cell(length(BFid),1);

% for curBF = 1:length(BFid)
%     ix = find(strcmp(BFlab,BFid(it)));
%     avgFR{curBF} = zeros(1,length(spkDat.fr(1,:)));
%     for curMW = 1:length(ix)
%         avgFR{curSesh}{curBF} = avgFR{curSesh}{curBF} + spkDat.fr(ix(curMW),:);
%     end;
% end;

return;