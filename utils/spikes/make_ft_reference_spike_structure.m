%%
restoredefaultpath;
addpath('~rouxf/fieldtrip-20161009/');
ft_defaults;

cd ~rouxf/fieldtrip-20161009/private/;
%cd ~rouxf/fieldtrip-20161009/external/neuralynx_v3/;

cfg = [];
cfg.dataset             = '/media/rouxf/rds-share/iEEG_DATA/MICRO/P04/fvsPEM/2016-10-18_15-30-26/CSC_L14.ncs';	
cfg.output              = '~rouxf/test/test.nse';
cfg.dataformat          = 'neuralynx_nse';
cfg.method              =  'zthresh';
cfg.interactive         = 'no';
cfg.timestampdefinition = 'orig';
    
cfg = ft_spikedetection(cfg);

%%
fn = '~rouxf/test/test.nse';

[spike] = ft_read_spike(fn);

%%
spike.cfg.forbidden = [];

%%
 
cfg =[];
cfg.channel = 'all';
cfg.method          ='kmeans';
cfg.feedback        ='no';
%cfg.kmeans          = [];
cfg.kmeans.aantal = 5;

[spike] = ft_spikesorting(cfg, spike)
