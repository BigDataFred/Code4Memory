function [idx] = getNlxHeaderParam(hdr,param)

x = regexp(hdr,['-',param]);
k = 0;
for it = 1:length(x)
    
    if ~isempty(x{it})
        k = k+1;
        idx(k) = it;
    end;
    
end;