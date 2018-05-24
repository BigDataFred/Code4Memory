[m,sd,r] = computeWVFstats(x)

%%
tmp = mean(x,2);
r = zeros(size(x,2),1);
for it = 1:size(x,2)
    r(it) = corr(x(:,it),tmp,'Type','Spearman');
end;

[pix] = find(sign(r)==1);

ppct = length(pix)/length(r);
m=tmp;
sd = std(x,0,2);