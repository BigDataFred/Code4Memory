function FFT_GUI()

restoredefaultpath;
addpath('/media/rouxf/rds-share/Common/fieldtrip-20170115/');
ft_defaults;

%%

[filen,pathn] = uigetfile({'*.*','All Files'});

%%
cfg                     = [];
cfg.dataset             = [pathn,filen];
%cfg.channel             = {'L Ent_4'};
cfg.trl                 = [1 1024*45 0];

[data] = ft_preprocessing(cfg);

%%
cfg                     = [];
cfg.lpfilter            = 'yes';
cfg.lpfreq              = data.fsample/2;
cfg.lpfilttype          = 'fir';

[data] = ft_preprocessing( cfg , data );
%%
cfg                     = [];
cfg.resamplefs          = 250;
cfg.detrend             = 'no';

[data] = ft_resampledata( cfg , data );

%%

cfg = [];
cfg.detrend             = 'yes';
cfg.demean              = 'yes';

[data] = ft_preprocessing( cfg , data );

%%
cfg                     = [];
cfg.method              = 'mtmfft';
cfg.taper               = 'dpss';
cfg.tapsmofrq           = 1;
cfg.pad                 = 'nextpow2';

[pow] = ft_freqanalysis( cfg , data );

%%
nch = length(pow.label);
figure;
for it = 1:nch
    subplot(ceil(nch/10),10,it);
    plot(pow.freq,pow.powspctrm(it,:));
    axis tight;
end;