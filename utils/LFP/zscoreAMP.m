function [zDat] = zscoreAMP(dat,selIdx)

zDat = cell(1,length(dat));
for it = 1:length(dat)
    
    x = dat{it};
    if nargin <2 || isempty(selIdx)
        M = repmat(mean(x,1),[size(x,1) 1]);
        SD = repmat(std(x,0,1),[size(x,1) 1]);
        z = (x-M)./SD;
        zDat{it} = z;
    else
        M = repmat(mean(x(selIdx,:),1),[size(x,1) 1]);
        SD = repmat(std(x(selIdx),0,1),[size(x,1) 1]);
        z = (x-M)./SD;
        zDat{it} = z;
    end;
end;