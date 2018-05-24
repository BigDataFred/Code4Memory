function [corr_S] = correct_one_over_f(S, f)

log_f = log10(f);
log_S = 10*log10(S);

dim = size(log_f);
if dim(2) > dim(1)
    log_f = log_f';
end;

dim = size(log_S);
if dim(2) > dim(1)
    log_S = log_S';
end;

%[b] = regress(log_S,[ones(size(log_S)) log_f]);
[b] = robustfit(log_f,log_S);
yp = b(1)+b(2)*log_f;

delta = log_S - yp;

%zsc = (delta-mean(delta))./std(delta);
%ix = find(abs(zsc)>1.5);

%x = 1:length(delta);
%y = yp;
%x(ix) = [];
%y(ix) = [];

%[yi] = interp1(x,y,ix);

%yp(ix) = yi;
%yp(isnan(yp)) = 0;

%delta = log_S - yp;

corr_S = 10.^delta;

return;
