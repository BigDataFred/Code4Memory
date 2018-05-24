function mergedCluster = merge_clusters(mergedCluster,clusterDat,clusterIdx)

S = fieldnames(clusterDat);
   
if isempty(mergedCluster)
    
    for it = 1:length(S)
        if iscell(clusterDat.(S{it}))
            dum = clusterDat.(S{it})(clusterIdx);
            clusterDat.(S{it}) = {};
            x = [];
            for jt = 1:length(dum)
                x = [x squeeze(dum{jt})];
            end;
            clusterDat.(S{it}) = x;
        end;
    end;
    mergedCluster = clusterDat;
    
else
    
    for it = 1:length(S)
        if iscell(mergedCluster.(S{it}))            
            dum =  {mergedCluster.(S{it}) clusterDat.(S{it})(clusterIdx)};
            mergedCluster.(S{it}) = {};
            x = [];
            for jt = 1:length(dum)
                x = [x squeeze(dum{jt})];
            end;
            mergedCluster.(S{it}) = x;
        end;
    end;
    
end;


return;