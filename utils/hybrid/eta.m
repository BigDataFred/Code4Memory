function [eta] = eta(dt,esmp,lfp_data,params)    
if nargin ==0
    dt = 0.006;
end;

c1 = 0;
eta = zeros(dt*params.Fs+1,length(esmp));
for it = 1:length(esmp)
    if sign(esmp(it)-dt/2*params.Fs)==1 && (esmp(it)+dt/2*params.Fs <=length(lfp_data))
        c1 =c1+1;
        eta(:,c1) = lfp_data(floor((esmp(it)-dt/2*params.Fs:esmp(it)+dt/2*params.Fs)+eps));
    end;
end;
eta(:,c1+1:end) = [];