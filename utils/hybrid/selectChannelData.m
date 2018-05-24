function [outDat] = selectChannelData(selIdx, varargin)

x = varargin;
x = x{:};
outDat = cell(1,length(x));
for it = 1:length(x)
    if size(x{it},2) > size(x{it},1)
        outDat{it} = x{it}(:,selIdx);
    elseif size(x{it},2) == size(x{it},1)
        outDat{it} = x{it}(selIdx,selIdx);
    else
        outDat{it} = x{it}(selIdx,:);
    end;
end;