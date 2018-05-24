function [TFR] = time_freq_EM( LFP_data, varargin )

TFR = {};
%%
restoredefaultpath;
addpath(['~rouxf',filesep,'fieldtrip-20161009',filesep]);
ft_defaults;

%% remove DC, detrend & 1/F
cfg                 = [];
% cfg.detrend         = 'yes';
% cfg.demean          = 'yes';
cfg.derivative      = 'yes';

[LFP_data] = ft_preprocessing(cfg, LFP_data );

%% Lower frequencies
if strcmp(varargin{1},'low') || strcmp(varargin{1},'low&high')
        %%
        cfg                 = [];
        if ~isempty(varargin{2})
            cfg.channel = LFP_data.label(varargin{2});
        end;
        cfg.method          = 'mtmconvol';
        cfg.taper           ='hanning';
        cfg.tapsmofrq       = 1;
        cfg.pad             = 'maxperlen';
        cfg.foi             = 2:25;
        cfg.t_ftimwin       = ones(1,length(cfg.foi));
        cfg.toi             = LFP_data.time{1}(1):0.25:LFP_data.time{1}(end);
        cfg.keeptrials      = 'yes';
        
        [TFRL] = ft_freqanalysis(cfg , LFP_data );
        
        cfg2                 = [];
        cfg2.latency         = [TFRL.time(1)+cfg.t_ftimwin(1) TFRL.time(end)-cfg.t_ftimwin(1)];
        
        [TFRL]              = ft_selectdata(cfg2, TFRL );
        
        TFR{end+1} = TFRL;
end;

%% Higer frequencies
if strcmp(varargin{1},'high') || strcmp(varargin{1},'low&high')        
        %%
        cfg                 = [];
        if ~isempty(varargin{2})
            cfg.channel = LFP_data.label(varargin{2});
        end;
        cfg.method          = 'mtmconvol';
        cfg.taper           ='dpss';
        cfg.tapsmofrq       = 15;
        cfg.pad             = 'maxperlen';
        cfg.foi             = 25:130;
        cfg.t_ftimwin       = 0.2*ones(1,length(cfg.foi));
        cfg.toi             = LFP_data.time{1}(1):0.05:LFP_data.time{1}(end);
        cfg.keeptrials      = 'yes';
        
        [TFRH] = ft_freqanalysis(cfg , LFP_data );
        
        cfg2                 = [];
        cfg2.latency         = [TFRH.time(1)+cfg.t_ftimwin(1) TFRH.time(end)-cfg.t_ftimwin(1)];
        
        [TFRH]              = ft_selectdata(cfg2, TFRH );
        
        TFR{end+1} = TFRH;
end;

return;