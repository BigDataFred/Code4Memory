%%
ts = spike1.timestamp{1}-  params.hdr.FirstTimeStamp;

ttl_idx = [params.event(:).value];
ttl_idx = find(ttl_idx ==7);
picon = [params.event(:).timestamp];
picon = picon(ttl_idx) -  uint64(params.hdr.FirstTimeStamp);

pic_sel = find(LogDat.stimID == LogDat.ID(27));

figure;
hold on;
for it = 1:length(pic_sel)
    
    x = ts - int64(picon(pic_sel(it)));
    x = double(x).*1e-6;
    
    x = x(find(x>=-params.pre & x <=params.post));
    
    x = [x;x];
    y = it*ones(1,length(x));
    y = [y-.5;y+.5];
    
    h = [];
    h = line(x,y);
    set(h,'Color','k');
    
end;