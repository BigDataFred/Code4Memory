%%
addpath('/home/rouxf/tbx/fieldtrip-20170618/');
ft_defaults;
addpath('/home/rouxf/prj/Bham/code/mcode/utils/LFP/');

%%
p2d = '/media/rouxf/rds-share/iEEG_DATA/MICRO/P02/fVSpEM/2016-07-13_14-14-25/lfp_dat/';
files = dir([p2d,'*_stimlocked.mat']);

lfpDat = cell(1,length(files));
for it = 1:length(files)
    fprintf([num2str(it),'/',num2str(length(files))]);
    load([p2d,files(it).name]);
    [lfpDat{it}] = save_data{1}{1}{1};
    fprintf('\n');
end;

[lfpDat] = ft_appenddata([],lfpDat{:});
lfpDat.fsample = 1/(lfpDat.time{1}(2)-lfpDat.time{1}(1));

%%
BFlabel = {}; 
for it = 1:length(lfpDat.label);
    BFlabel(it) = {lfpDat.label{it}(1:regexp(lfpDat.label{it},'\d')-1)};
end

mwID = unique(BFlabel);

mwAVG = {};
for it = 1:length(mwID)
    
    selIx = find(strcmp(BFlabel,mwID{it}));
    
    cfg                     = [];
    cfg.channel             = lfpDat.label(selIx);
    cfg.avgoverchan         = 'yes';
    
    [mwAVG{it}] = ft_selectdata( cfg , lfpDat );
    
end;

[mwAVG] = ft_appenddata( cfg , mwAVG{:} );

%%
cfg = [];
cfg.viewmode = 'vertical';

ft_databrowser(cfg,mwAVG);

%%
cfg = [];
cfg.hpfilter = 'yes';
cfg.hpfreq = [15];

[filtDat] = ft_preprocessing(cfg,mwAVG);

for it = 1:length(filtDat.label);filtDat.label(it) = {['dum_',num2str(it)]};end;

%%
filtDat = ft_appenddata([], mwAVG, filtDat );

%%
cfg = [];
cfg.viewmode = 'vertical';

ft_databrowser(cfg,filtDat);

%%
cfg                     = [];
cfg.channel             = {'dum*'};

[filtDat] = ft_selectdata( cfg , filtDat );

%%
rejTrl = [];
for it = 1:length(lfpDat.trial)
    
    zsc = zscore(filtDat.trial{it}')';
    zsc = zsc >8;
    zsc = sum(zsc,1);
    if any(zsc >= 1)
        rejTrl = [rejTrl it];
    end;
    
end;

%%
cfg                     = [];
cfg.trials              = rejTrl;

[rejDat] = ft_selectdata( cfg , mwAVG );

%%
cfg                     = [];
cfg.trials              = setdiff(1:length(mwAVG.trial), rejTrl);

[cleanDat] = ft_selectdata( cfg , mwAVG );

%%
cfg                     = [];
cfg.viewmode            = 'vertical';

ft_databrowser(cfg,cleanDat);

%%
cfg                     = [];
cfg.viewmode            = 'vertical';

ft_databrowser(cfg,rejDat);

%%
cfg                     = [];
cfg.method              = 'wavelet';
cfg.width               = 7;
cfg.gwidth              = 5;
cfg.toi                 = lfpDat.time{1}(1):0.025:lfpDat.time{1}(end);
cfg.foi                 = 0:300;
cfg.keeptrials          ='yes';
cfg.output              = 'pow';

[pow] = ft_freqanalysis( cfg , lfpDat );