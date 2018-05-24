function [Y] = sortedEVD(X)

[nchan,nsmp,ntrl] = size(X);

Y = zeros(nchan,nsmp,ntrl);
pct = zeros(nchan,ntrl);
for jt = 1:ntrl
    
    x = squeeze(X(:,:,jt));
    m = mean(x,2)*ones(1,size(x,2));
    x = x-m;
    x = x';
    R = x'*x;
    
    [v,d] = eig(R);
        
    v(cumsum((diag(d)./sum(diag(d))*100))>=2,:) = 0;
    
    y = (v'*x')';
    Y(:,:,jt) = y';
end;

return;