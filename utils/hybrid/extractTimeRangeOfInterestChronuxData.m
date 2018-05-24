function [dataOut] = extractTimeRangeOfInterestChronuxData(selIdx,toi,dataIn)

pre = toi(1);
post = toi(2);
dataOut = cell(1,length(dataIn));
if isstruct(dataIn{1})
    for it = 1:length(dataIn)
        
        x = dataIn{it};
        for jt = 1:length(x)
            x(jt).times = x(jt).times(x(jt).times>=pre & x(jt).times <=post);
        end;
        dataOut{it} = x;
        
    end;
else    
    parfor it = 1:length(dataIn)
        
        dataOut{it} = dataIn{it}(selIdx,:);
        
    end;        
end;