function [Fval,Critval,pval] = computeFtest( design , dat)

condID= unique(design);
ncond = length(condID);
nrepl=0;
for condindx=1:ncond
    nrepl=nrepl+length(find(design==condID(condindx)));
end;
if nrepl<length(design)
  error('Invalid specification of the independent variable in the design array.');
end;
if nrepl<=ncond
    error('The must be more trials/subjects than levels of the independent variable.');
end;
dfnum = ncond - 1;
dfdenom = nrepl - ncond;

nsmpls = size(dat,1);

% compute the statistic
nobspercell = zeros(1,ncond);
avgs = zeros(1,ncond);
pooledvar = zeros(1,1);
for condindx=1:ncond
    sel=find(design==condindx);
    nobspercell(condindx)=length(sel);
    avgs(condindx)=mean(dat(sel));
    pooledvar = pooledvar + nobspercell(condindx)*var(dat(sel),1);
end;
pooledvar = pooledvar/dfdenom;
globalavg = mean(mean(dat));
mseffect = ((avgs-repmat(globalavg,1,ncond)).^2)*nobspercell'./dfnum;

Fval = mseffect./pooledvar;
Critval = finv(1-0.05,dfnum,dfdenom);
pval = 1-fcdf(Fval,dfnum,dfdenom);