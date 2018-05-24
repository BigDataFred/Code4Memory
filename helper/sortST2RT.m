function [dat] = sortST2RT(RT,dat)
%%
[~,s_idx] = sort(RT);

x1 = []; x2 = []; x3 = [];
for kt = 1:length(s_idx)
    
    
    ix = find(dat.trial{1} == s_idx(kt));
    
    x1 = [x1 dat.trial{1}(ix)];
    x2 = [x2 dat.time{1}(ix)];
    x3 = [x3 kt*ones(1,length(ix))];
end;

if length(dat.trial{1}) ~= length(x1) || length(dat.time{1}) ~= length(x2)
    error('number of elements is out of range');
end;

dat.trial{1} = x1;
dat.trial2{1} = x3;
dat.time{1} = x2;
