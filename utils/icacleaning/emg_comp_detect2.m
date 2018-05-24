function [emg_idx] = emg_comp_detect2(comp)

%% high-pass filter and hilbert transform
cfg = [];
cfg.hpfilter = 'yes';
cfg.hpfreq = 30;
cfg.hilbert = 'abs';

[comp] = ft_preprocessing(cfg,comp);
%%
nsamp = zeros(length(comp.trial),1);
for it = 1:length(nsamp)
    nsamp(it) = length(comp.trial{it});
end;

concat = zeros(length(comp.label),sum(nsamp));
idx = 1:length(comp.trial{1});
for it = 1:length(nsamp)
    concat(:,idx) = comp.trial{it};
    if sign(it-length(comp.trial))==-1
        idx = idx(end)+1:idx(end)+length(comp.trial{it+1});
    end;
end;
%% get the stationarity acrros the data
win = 1:comp.fsample*1;
stp = length(win)*0.75;
nsty = zeros(length(comp.label),round(size(concat,2)/stp));
for jt = 1:length(comp.label)
    k=0;
    win = 1:comp.fsample*1;
    stp = length(win)*0.75;
    for it = 1:round(size(concat,2)/stp)
        
        if sign(size(concat,2)-win(end))==1
            k=k+1;
            nsty(jt,k) = std(concat(jt,win),[],2);        
            win = win+stp;
        end;
    end;

end;
nsty(:,k+1:end)=[];
%%
dum = struct;
dum.time{1} = 1:size(nsty,2);
dum.fsample = 1/(dum.time{1}(2)-dum.time{1}(1));
dum.trial{1} = nsty;
dum.label = comp.label;

cfg = [];
cfg.method = 'mtmfft';
cfg.pad = 'maxperlen';
cfg.output = 'pow';
cfg.taper = 'dpss';
cfg.tapsmofrq = .1;

pow = ft_freqanalysis(cfg,dum);

dum = struct;
dum.time{1} =pow.freq;
dum.fsample = 1/(dum.time{1}(2)-dum.time{1}(1));
dum.trial{1} = (pow.powspctrm);
dum.label = comp.label;


cfg = [];
cfg.method = 'mtmfft';
cfg.pad = 'maxperlen';
cfg.output = 'pow';
cfg.taper = 'dpss';
cfg.tapsmofrq = 2;

pow = ft_freqanalysis(cfg,dum);

df=sum(pow.powspctrm,2);
df = df-mean(df);
df = df./std(df);

[df,s_idx] = sort(df);

[emg_idx] = s_idx(find(sign(df-0.1)==1));


