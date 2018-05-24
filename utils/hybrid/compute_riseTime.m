function [rt] = compute_riseTime(wvf)

for it =1:size(wvf,1)
    
    x = abs(wvf(it,:));
    m = max(x);
    x = x./m;
    
    ix2 = find(x == max(x));
    ix1 = find(x(ix2:-1:1)<0.2);
    ix3 = find(x(ix2:end)<0.2);
    
    ix1 = ix1(1);
    ix3 = ix2+ix3(1);
    
    ix = [ix1:ix3];
    
    rt(it) = length(ix);
    
   
    
end;
