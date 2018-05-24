function [trl] = computeTRLmatrix(pre,post,events,lfpTime)    
if nargin ==0
    pre = 5;
    post = 7;
end;

[trl] = zeros(size(events,1),3);
for it = 1:size(events,1)
    ix = find( lfpTime >= (events(it,1)-pre) & lfpTime <= (events(it,1)+post) );
    [trl(it,:)] = [min(ix) max(ix) -pre*Fs];
end;