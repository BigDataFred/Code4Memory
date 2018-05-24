function [ITC,ERPimg,tAx] = computeITC4LFP(LFPdat,foi,toi,chanLab)


cfg                     = [];
cfg.channel             = chanLab;

[LFPdat] = ft_selectdata( cfg , LFPdat);


cfg                     = [];
cfg.keeptrials          = 'yes';
cfg.preproc.lpfilter    = 'yes';
cfg.preproc.lpfreq      = [30]; 
cfg.preproc.lpfilttype  = 'but';
cfg.preproc.lpfiltord   = 4;


[ERPimg] = ft_timelockanalysis( cfg , LFPdat );

cfg                     = [];
cfg.latency             = [toi(1) toi(end)];


[ERPimg] = ft_selectdata( cfg , ERPimg );

tAx = ERPimg.time;
ERPimg = ERPimg.trial;

cfg                     = [];
cfg.channel             = chanLab;
cfg.method              = 'wavelet';
cfg.foi                 = foi;
cfg.toi                 = LFPdat.time{1}(1):0.025:LFPdat.time{1}(end);
cfg.width               = 5; 
cfg.output              = 'fourier';
cfg.pad                 = 'nextpow2';


[phi] = ft_freqanalysis( cfg , LFPdat );

cfg                     = [];
cfg.latency             = [toi(1) toi(end)];


[phi] = ft_selectdata( cfg , phi );

chck = ~isnan(angle(squeeze(phi.fourierspctrm(1,1,1,:))));

F = phi.fourierspctrm;
tmp = F./abs(F);
ITC = abs(mean(tmp,1));
ITC = squeeze(ITC);
