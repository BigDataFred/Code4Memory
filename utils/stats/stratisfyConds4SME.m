function [stratDat,design2] = stratisfyConds4SME(FR,design)
%%
origDat = FR;
condID = unique(design);

ncond=[];
for jt = 1:length(condID)
    ncond(jt) = length(find(design == condID(jt)));
end;
[~,ix] = min(ncond);

%%
design2 = [ones(1,ncond(ix)-1) 2*ones(1,ncond(ix)-1) 3*ones(1,ncond(ix)-1)];

%%
stratDat = {};
for jt = 1:length(condID)
    sel = find(design ==jt);
    sel = sel(randperm(length(sel)));
    sel = sel(1:ncond(ix)-1);
    stratDat{jt} = origDat(sel);
end;

