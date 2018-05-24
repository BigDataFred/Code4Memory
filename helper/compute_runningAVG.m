function [ra] = compute_runningAVG(data,w)
%%
n = length(data);
data = [zeros(1,w) data' zeros(1,w)];
c = 0;
ra = zeros(1,length(w+1:length(data)-w));
for it = w+1:length(data)-w
    c =c +1;
    x = data(it-w:it+w);
    ra(c) = mean(x(x~=0));
end;

if length(ra) ~=n
    error('vector length is out of range');
end;
return;
