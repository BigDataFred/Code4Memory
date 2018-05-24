%%
figure;
hold on;
id = unique(dumm.trial{1});%find(logpicident == find(strcmp(picnamev,'f_ppp243.jpg')));%
x = cell(1,length(id));
for it = 1:length(id)
    
    sel_idx =find(dumm.trial{1} == id(it));
    
    [dum] = double(dumm.timestamp{1}(sel_idx)-int64(event(ttl_idxc(it)).timestamp));
    dum = dum./1e6;
    
    dum = dum(find(dum >= -.5 & dum <=1));
    
    if abs(min(dum)) >0.5
        error('tata');
    end;
       
    if abs(max(dum)) >1
        error('toto');
    end;
    
    if isempty(dum)
        error('titi');
    end;
    
    x{it} = dum;
    
    y = repmat(it,[1 length(x{it})]);
    
    line([x{it};x{it}],[y-.5;y+.5],'Color',[0 0 0]);
    
end;
axis tight;axis xy;
xlim([-.5 1]);
%%
n1 = [];
n2 = [];
r = [];
for it = 1:length(x)
    
    n1(it) = length(x{it});
    n2(it) = length(chck.x{it});
    try
        r(it) = corr(x{it}',chck.x{it}');
    catch
        r(it) = NaN;
    end;
end;
%%
sel_idx =find(logpicident == find(strcmp(picnamev,'f_ppp243.jpg')));%
figure;
hold on;
for it = 1:length(sel_idx)
    d = x{sel_idx(it)};
    r2(it)= corr(d',chck.x{sel_idx(it)}');
    y = repmat(it,[1 length(d)]);
    
    line([d;d],[y-.5;y+.5],'Color',[0 0 0]);
end;