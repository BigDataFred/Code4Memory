function [data] = remove_cluster_field(data,sel_idx)

S = fieldnames(data);

for it = 1:length(S)
    if iscell(data.(S{it}))
        data.(S{it})(sel_idx) = [];
    end;
end;