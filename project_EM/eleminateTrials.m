function [outDat] = eleminateTrials(inDat,delIx)

for it = 1:length(inDat)

    inDat{it}(delIx) = [];

end;
outDat = inDat;

return;