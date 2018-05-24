%%

S = load('~rouxf/CSC_RP6_spikes.mat')

dumm = [];
dumm.label = {'CSC_RA6'};
dumm.timestamp{1} = uint64(S.index);
dumm.waveform{1}(1,:,:) = uint64(S.spikes(:,1:2:end)');
dumm.unit{1} = zeros(1,length(S.index));
dumm.hdr = hdr;
dumm. dimord = '{chan}_lead_time_spike';

dumm

cfg                         = [];
%cfg.trl                     = trl2(ix,:);
cfg.trl                     = trl;
%cfg.trlunit                 = 'timestamps';
%cfg.timestampspersecond     = hdr.Fs*hdr.TimeStampPerSample;
cfg.trlunit                 = 'samples';
cfg.hdr                     = hdr;

[dumm] = ft_spike_maketrials(cfg,spike3);

save('testdataMATT_CSC_RA6.mat','dumm','LogDat')
%%
 f = [];
f.p2logf = '/media/rouxf/rds-share/EXP/EMpairs_v4_2016-1007/log/Tuning/';
f.logf = Tlogfile{1};

p = [];
p.ncols = 8;
[LogDat] = getNewLogDataTune(f,p);
%%
buf = LogDat.dat(1:end,1);
logpicident =[];
for n = 1:length(buf)
    logpicident(n) = str2num(buf{n});
end;
numpics = length(unique(logpicident));

picnamev = [];
for n = 1:numpics
    f = find(logpicident == n,1,'first');
    picnamev{n} = LogDat.dat{f,2};
end;

titles = {};
for p = 1:length(picnamev)
    [a,titles{p},b] = fileparts(picnamev{p});
end;
%%

cfg = [];
cfg.latency             = [0.1 0.9];
cfg.trials = find(logpicident == find(strcmp(picnamev,'a_sea28.jpg')));
cfg.plotselection   =  'yes';

figure;
subplot(4,1,1:3)
ft_spike_plot_raster(cfg,dat{30}.spike3);
axis xy;

cfg                     = [];
cfg.binsize             = 0.011;
cfg.outputunit          = 'rate';
cfg.latency             = [0.1 0.9];
cfg.trials = find(logpicident == 27);

[psth] = ft_spike_psth(cfg,dat{30}.spike3);

subplot(4,1,4);
plot(psth.time,psth.avg);

%%
r =[];
s=[];
for p = 1:numpics
    
    f = find(logpicident == p);
    
    %%
    cfg                     = [];
    cfg.outputunit          = 'spikecount';
    cfg.keeptrials          = 'yes';
    cfg.binsize             = 0.011;
    cfg.trials              = [];
    cfg.trials              = f;
    cfg.latency             = [-.5 1];
    
    [rate]              = ft_spike_psth(cfg,dumm);
    
    for z = 1:size(rate.trial,1)        
        rate.trial(z,:,:) = conv(squeeze(rate.trial(z,:,:)),gausswin(5),'same');
    end;
    
    cfg                     = [];
    cfg.latency             = [-.35 -.05];
    cfg.avgovertime         = 'no';
    
    [rateb]                   = ft_selectdata(cfg,rate);
    
    rate.trial = rate.trial.*rate.fsample - median(mean(rateb.trial.*rate.fsample,3),1);
    rate.avg = squeeze(mean(rate.trial,1));
    
    cfg                     = [];
    cfg.latency             = [.1 .9];
    cfg.avgovertime         = 'no';
    
    [rate]                   = ft_selectdata(cfg,rate);
    
        
    
    r(p,:) = [mean(rate.avg) std(squeeze(mean(rate.trial,3)))/sqrt(size(rate.trial,1))];
    
    [~,h,stats] = signtest(squeeze(mean(rate.trial.*rate.fsample,3)),0);
    if h ==1 && sign(stats.sign)==1
        s(p) = 1;
    end;
    
end;