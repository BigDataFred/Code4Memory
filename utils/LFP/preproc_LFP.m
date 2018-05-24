function [CSC_preproc] = preproc_LFP(dataset,trl)


%% read the trial segments and low pass
cfg1 = [];
cfg1.continuous      = 'yes';
cfg1.detrend         = 'yes';
cfg1.demean          = 'yes';
cfg1.dataset         = dataset;
cfg1.trl             = trl;

%     cfg1.lpfilter        = 'yes';
%     cfg1.lpfreq          = 340;
%     cfg1.lpfilttype      = 'fir';

[dum]                = ft_preprocessing(cfg1);

%     %% downsample
%     cfg2 = [];
%     cfg2.resamplefs     = 1360;
%
%     [dum]               = ft_resampledata(cfg2,dum);

%% do preprocessing
cfg3                 = [];
%     cfg3.dftfilter       = 'yes';
%     cfg3.dftfreq         = [50:50:200];
%     cfg3.bpfilter        = 'yes';
%     cfg3.bpfreq          = [0.5 170];
%     cfg3.bpfilttype      = 'fir';
%     cfg3.padtype         = 'mirror';
%     cfg3.padding         = round(length(dum.time{1})/dum.fsample*3.5);

[dum]                = ft_preprocessing(cfg3,dum);

[CSC_preproc] = dum;
clear_par(dum);

