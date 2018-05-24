function [pvalPL,zPL,PL] = computeSpkPL(spkDat,Phi)

[N,M] = size(spkDat);

[dims] = size(Phi);% 1=rpt,2=freq,3=time
pval  = zeros(1,dims(2));
z     = zeros(1,dims(2));
pl     = zeros(1,dims(2));
nfreq = dims(2);
parfor it = 1:nfreq% loop over freqs
    fprintf([num2str(it),'/',num2str(nfreq)]);
    phi = squeeze(Phi(:,it,:,:))';%transpose to time x rpt
    [n,m] = size(phi);
    if n~=N || m~=M
        error('matrix dimensions must lign up');
    end;
    ix = find(spkDat);
    if length(ix) >= 50                
        p = phi(ix);
        pl(it) = abs(sum(exp(1i*p)))/length(p);
        [pval(it), z(it)] = circ_rtest(p);
    end;
    fprintf('\n');
end;

pvalPL = pval;
zPL = z;
PL = pl;