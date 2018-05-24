function [lfp_dat] = load_lfp_data(lfp_files,path_2_lfps,trl_sel)

dat = cell( 1 , length( lfp_files ) );
parfor it = 1:length(dat)
    
    dum = load([path_2_lfps,lfp_files(it).name]);%loads the data for each LFP channel
    dat{it} = dum.save_data{1}{1}{1};
    dat{it}
end;

cfg                     = [];

[lfp_dat] = ft_appenddata( cfg , dat{:} );%concatenate all the LFP channels together
clear dat;

% select the trials and the time range of interest
cfg                     = [];
cfg.trials              = trl_sel;
cfg.latency             = [-2 4];

[lfp_dat] = ft_selectdata( cfg, lfp_dat);


cfg                     = [];
cfg.detrend             = 'no';
cfg.resamplefs          = 1000;

[lfp_dat] = ft_resampledata( cfg , lfp_dat );

% % put data into a conveniant format
% cfg                     = [];
% cfg.keeptrials          = 'yes';
% 
% [lfp_dat] = ft_timelockanalysis( cfg, lfp_dat);