%%
addpath('/media/rouxf/rds-share/Common/fieldtrip-20170108/');
ft_defaults;

%%
p2f = '/media/rouxf/rds-share/iEEG_DATA/MICRO/P04/fvSpEM/';
chck2 = dir(p2f);
chck2(1:2) = [];

%%
pow1 = cell(1,length(chck2));
pow2 = cell(1,length(chck2));

for zt = 1:length( chck2 )
    
    p2f2 = [p2f,chck2(zt).name,'/log_dat/'];
    
    fn = dir([p2f2,'*.mat']);
    load( [p2f2,fn.name] );
    
    sel = [];sel(1,:) = LogDat2.idx(1,:);for it = 2:size(LogDat2.idx,1);sel(it,:) = LogDat2.idx(it,:)+LogDat2.idx(it-1,2);end;
    sel2 = [];
    for it = 1:size(sel,1)
        sel2 = [sel2 sel(it,1):sel(it,2)];
    end;
    
    %%
    p2f2 = [p2f,chck2(zt).name,'/lfp_dat/'];
    files = dir([p2f2,'*.mat']);
    
    
    chck = regexp({files(:).name},'CSC_L');
    nL = [];
    ixL = [];
    k = 0;
    for it = 1:length(chck)
        
        if ~isempty(chck{it})
            k = k+1;
            ix2 = regexp(files(it).name(chck{it}:end),'_');
            ix2 = ix2(2);
            ixL(k) = it;
            nL(k) = str2double(files(it).name(chck{it}+5:chck{it}+ix2-2));
        end;
    end;
    
    chck = regexp({files(:).name},'CSC_R');
    nR = [];
    ixR = [];
    k = 0;
    for it = 1:length(chck)
        
        if ~isempty(chck{it})
            k = k+1;
            ix2 = regexp(files(it).name(chck{it}:end),'_');
            ix2 = ix2(2);
            ixR(k) = it;
            nR(k) = str2double(files(it).name(chck{it}+5:chck{it}+ix2-2));
        end;
    end;
    
    [~,s_idx] = sort(nL);
    ixL = ixL(s_idx);
    
    [~,s_idx] = sort(nR);
    ixR = ixR(s_idx);
    
    files = files([ixL ixR]);
    
    %%
    dat = cell( 1 , length(files) );
    for it = 1:length(files)
        
        load([p2f2,files(it).name]);
        
        dat{it} = save_data{1}{1}{1};
        clear save_data;
        
    end;
    
    %%
    dat                 = ft_appenddata( [] , dat{:} );
    
    %%
    cfg                 = [];
    cfg.trials          = sel2;
    
    [dat]               = ft_selectdata( cfg, dat );
    
    %%
    cfg                 = [];
    cfg.lpfilter        = 'yes';
    cfg.lpfreq          = 30;
    
    [dum]               = ft_preprocessing( cfg , dat );
    
    %%
    cfg                 = [];
    cfg.latency         = [-.5 3];
    
    [dum]               = ft_selectdata( cfg, dum );
    
    %%
    cfg                 = [];
    
    [tlck] = ft_timelockanalysis( cfg , dum );
    
    %%
    cfg                 = [];
    cfg.demean          = 'yes';
    cfg.detrend         = 'yes';
    
    [dum]               = ft_preprocessing( cfg , dat );
    
    cfg                 = [];
    cfg.derivative      = 'yes';
    
    [dum]               = ft_preprocessing( cfg , dum );
    
    cfg                 = [];
    cfg.latency         = [-1 0];
    
    [dum]               = ft_selectdata( cfg, dum );
    
    %%
    pad = zeros(length(dum.label),length(dum.time{1})*1);
    
    for it = 1:length(dum.trial)
        
        dum.trial{it} = [pad dum.trial{it} pad];
        dum.time{it}  = 0:1/dum.fsample:((length(dum.trial{it})-1)/dum.fsample);
        
    end;
    
    %%
    n = length(dum.time{1})/dum.fsample;
    
    cfg                 = [];
    cfg.method          = 'mtmfft';
    cfg.pad             = 'maxperlen';
    cfg.taper           = 'dpss';
    cfg.tapsmofrq       = 1/n;
    cfg.output          = 'pow';
    cfg.keeptrials      = 'no';
    
    [pow]               = ft_freqanalysis( cfg , dum );
    
    
    cfg                 = [];
    cfg.frequency       = [0 30];
    
    [pow1{zt}]               = ft_selectdata( cfg , pow );
    
    %%
    cfg                 = [];
    cfg.demean          = 'yes';
    cfg.detrend         = 'yes';
    
    [dum]               = ft_preprocessing( cfg , dat );
    
    cfg                 = [];
    cfg.derivative      = 'yes';
    
    [dum]               = ft_preprocessing( cfg , dum );
    
    cfg                 = [];
    cfg.latency         = [0 3];
    
    [dum]               = ft_selectdata( cfg, dum );
    
    %%
    % pad = zeros(length(dum.label),length(dum.time{1})*5);
    %
    % for it = 1:length(dum.trial)
    %
    %     dum.trial{it} = [pad dum.trial{it} pad];
    %     dum.time{it}  = 0:1/dum.fsample:((length(dum.trial{it})-1)/dum.fsample);
    %
    % end;
    
    %%
    n = length(dum.time{1})/dum.fsample;
    
    cfg                 = [];
    cfg.method          = 'mtmfft';
    cfg.pad             = 'maxperlen';
    cfg.taper           = 'dpss';
    cfg.tapsmofrq       = 1/n;
    cfg.output          = 'pow';
    cfg.keeptrials      = 'no';
    
    [pow]               = ft_freqanalysis( cfg , dum );
    
    
    cfg                 = [];
    cfg.frequency       = [0 30];
    
    [pow2{zt}]               = ft_selectdata( cfg , pow );
end;

%%

cfg                     = [];
cfg.parameter           = 'powspctrm';

[pow1]                    = ft_appendfreq(cfg,pow1{:});
[pow2]                    = ft_appendfreq(cfg,pow2{:});

%%
a = [];
lm = [];

figure;
for it = 1:length(pow.label)
    subplot(6,8,it);
    a(it) = gca;
    hold on;
    plot(pow1.freq,squeeze(mean(pow1.powspctrm(:,it,:),1)),'r');
    plot(pow2.freq,squeeze(mean(pow2.powspctrm(:,it,:),1)));
    axis tight;
    lm = [lm min(mean(pow2.powspctrm(:,it,:),1)) max(mean(pow2.powspctrm(:,it,:),1))];
    title(pow2.label(it));
end;
set(a,'YLim',[min(lm) max(lm)]);