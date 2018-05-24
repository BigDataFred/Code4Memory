function [lfpDat] = sigProc4TimeFreqAnalysis(lfpDat)

movingwin = [1 .5];

for curMW = 1:length(lfpDat.chanLab)
    fprintf([num2str(curMW),'/',num2str(length(lfpDat.chanLab))]);
    lfpDat.LFPseg{curMW} = lfpDat.LFPseg{curMW}-ones(size(lfpDat.LFPseg{curMW},1),1)*mean(lfpDat.LFPseg{curMW},1);
    X = lfpDat.LFPseg{curMW};
%     parfor curTrl = 1:size(X,2)
%         x = X(:,curTrl)';
%         %f = filtfilt(x,1,[fliplr(x) x fliplr(x)]);
%         x = f(length(x)+1:length(x)*2);
%         x = locdetrend(x,lfpDat.Fs,movingwin);
%         X(:,curTrl)= x;
%     end;
    lfpDat.LFPseg{curMW} = X;
    fprintf('\n');
end;