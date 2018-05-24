
function [avgSxx1,avgSxx2,BFid] = avgLFPPowData(LFPpowDat,BFlab)

[BFid] = unique( BFlab );

[avgSxx1] = cell(length( BFid ),1);
[avgSxx2] = cell(length( BFid ),1);

for curBF = 1:length( BFid )
    [avgSxx1{curBF}] = zeros(size(LFPpowDat{1}.Sxx1));
    [avgSxx2{curBF}] = zeros(size(LFPpowDat{1}.Sxx2));
    MWix = find(strcmp(BFlab,BFid(curBF)));
    for curMW = 1:length( MWix )
        fprintf([num2str(curMW),'/',num2str(length( MWix ))]);
        LFPpowDat{MWix(curMW)}.Sxx1 = normalize(LFPpowDat{MWix(curMW)}.Sxx1);
        [avgSxx1{curBF}] = avgSxx1{curBF} + LFPpowDat{MWix(curMW)}.Sxx1;
        LFPpowDat{MWix(curMW)}.Sxx2 = normalize(LFPpowDat{MWix(curMW)}.Sxx2);
        [avgSxx2{curBF}] = avgSxx2{curBF} + LFPpowDat{MWix(curMW)}.Sxx2;
        fprintf('\n');
    end;
    [avgSxx1{curBF}] = avgSxx1{curBF}./length( LFPpowDat );
    [avgSxx2{curBF}] = avgSxx2{curBF}./length( LFPpowDat );
end;


