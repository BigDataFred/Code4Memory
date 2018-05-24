function [ttl_idx] = check_ttl_vs_log(log_params,ttl,ntrl)
%%
tc = zeros(length(log_params.idx),1);
for it = 1:length(log_params.idx)
    
    tc(it) = log_params.tc(log_params.idx(log_params.idx(it)));
    
end;

c = zeros(length(log_params.tc),1);
idx2 = zeros(ntrl,1);
for it = 1:ntrl
    
    idx = find(ttl == tc(it));
    
    c((tc(it)==log_params.tc)) = c((tc(it)==log_params.tc)) + 1;
    
    idx2(it) = idx(c((tc(it)==log_params.tc)));
    
end;
ttl_idx = idx2;

dx = diff([tc(1:ntrl) ttl(ttl_idx)'],[],2);
if sum(dx) ~=0
    error('trigger assignment must match');
end;