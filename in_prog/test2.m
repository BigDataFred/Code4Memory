%%
picev = log{46}.LogDat2.dat(:,2);
picid = unique(picev);

sel_idx = cell(1,length(picid));
for xt = 1:length(picid)
    sel_idx{xt} = find(strcmp(picev,picid(xt)));
end;
sel_idx = sel_idx';
%%
cfg = [];
cfg.latency = [.1 .9];
cfg.trials = find(log{30}.LogDat2.stimID == log{30}.LogDat2.ID(68));   

[dum] =  ft_selectdata(cfg,dat{30}.spike2);
%dum.trialtime = repmat([min(dum.time{1}) max(dum.time{1})],[6 1]);

figure;
ft_spike_plot_raster([],dum);
%%
x = dat{30}.spike2;

trl = sort(unique(x.trial{1}));
figure;
hold on;
for it = 1:length(trl)
    
    if ismember(it,1:size(x.trialtime,1))
        ix = find(ismember(x.trial{1},trl(it)));
        g = find(trl(it) == 1:size(x.trialtime,1));
        if g~=0
            plot(x.time{1}(ix),g*ones(1,length(ix)),'b.');    
        end;
    end;
end;
xlim([-.5 1]);
ylim([1 size(x.trialtime,1)]);