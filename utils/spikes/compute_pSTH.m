function [pSTH] = compute_pSTH(event,stts,dt,bw)

% Create a pool using default settings and disable the use of SPMD.
%parpool('SpmdEnabled', false);

dt = abs(dt);
if length(dt) ==1
    dt = [dt dt];
end;

nt = length(0:bw:sum(dt));
    
pSTH = zeros(length(event),nt);
for it = 1:length(event)
    
    x = event(it);
    xb = x-dt(1):bw:x+dt(2);
    
    n = zeros(length(xb),1);
    for jt = 1:length(xb)
        if jt < length(xb)
            n(jt) = length(find(stts >= xb(jt) & stts < xb(jt+1)));               
        end;
    end;
    pSTH(it,:) = n;
end;

%delete(gcp);

return