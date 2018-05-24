function [erpDat,TFR] = compute_ERP_and_TFR_macroData_EM(macroDat)

%% downsample the original data
cfg                     = [];
cfg.resamplefs          = 600;
cfg.detrend             = 'no';

[macroDat] = ft_resampledata( cfg , macroDat );

cfg                     = [];
cfg.lpfilter            = 'yes';
cfg.lpfreq              = 10;

[lpFilt] = ft_preprocessing( cfg , macroDat );

%% compute the ERP 4 encoding
cfg                     = [];

[dum] = ft_selectdata( cfg , lpFilt );

% cfg                     =[];
% cfg.method              = 'runica';
% 
% [dum] = ft_componentanalysis( cfg , dum);

cfg                     = [];
cfg.keeptrials          = 'yes';

[erpDat] = ft_timelockanalysis( cfg , dum );
clear dum;

%% compute TFR 4 encoding
cfg                     = [];
% cfg.detrend             = 'yes';
% cfg.demean              = 'yes';
%cfg.derivative          = 'yes';

[dum] = ft_preprocessing( cfg , macroDat );

cfg                     = [];
cfg.method              = 'mtmconvol';
cfg.pad                 = 'maxperlen';
cfg.taper               = 'dpss';
cfg.foi                 = [1:19 20:4:170];
cfg.toi                 = dum.time{1}(1):0.01:dum.time{1}(end);
cfg.keeptrials          = 'yes';

cfg.tapsmofrq(find(cfg.foi<20))     = 1;
cfg.tapsmofrq(find(cfg.foi>=20))    = 10;

cfg.t_ftimwin(cfg.foi<20)           = ones(1,length(find(cfg.foi<20)));
cfg.t_ftimwin(cfg.foi>=20)          = 0.25*ones(1,length(find(cfg.foi>=20)));

[TFR] = ft_freqanalysis( cfg , dum );
clear dum;

