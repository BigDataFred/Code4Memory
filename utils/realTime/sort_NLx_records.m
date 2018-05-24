function [ts,samp] = sort_NLx_records(ts,samp)

[v,ix] = sort(ts);

ts = ts(ix);
samp = samp(:,ix);

return;