function [AC]= compute_spiketrain_autocorr(st,varargin)
%%
if isempty( varargin )
    bw = .001;
    l = 1;
else
    bw = varargin{1};
    l  = varargin{2};
end;
% Create a pool using default settings and disable the use of SPMD.
%parpool('SpmdEnabled', false);

if size(st,2) > size(st,1)
    st = st';
end;

x = cell(length(st),1);
parfor it = 1:length(st)
    
    b = [st(it)-l st(it)+l];
    
    idx = find(st >= min(b) & st <= max(b));

    dt = diff([st(it)*ones(length(idx),1) st(idx)],[],2);
    
    b = linspace(min(b)-st(it),max(b)-st(it),2*1/bw+1);
    
    %x{it} = zeros(length(dt),length(b));
    x{it} = zeros(1,length(b));
    %for kt = 1:length(dt)
        for nt = 1:length(b)
            if nt < length(b)
                %x{it}(kt,nt) = length( find( dt-dt(kt) >= b(nt) & dt-dt(kt) < b(nt+1) ) );
                x{it}(nt) = length( find( dt >= b(nt) & dt < b(nt+1) ) );
            end;
        end;
    %end;
    x{it} = x{it}./(length(dt)*bw);
    
end;

b = linspace(-l,l,2*1/bw+1);

AC.trl = x;
X = zeros(1,length(b));
parfor it = 1:length(x)
    X = X+(sum(x{it},1));
end;
X = X./length(x);

AC.avg = X;

t = cell(1,length(x));
parfor it = 1:length(x)
    t{it} = b;
end;
AC.time = t;

AC.units = 's';
AC.bw = bw;

delete(gcp);

return;