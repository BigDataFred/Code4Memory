%%
addpath('/media/samba_share/RDS/Fred/code/mcode/toolboxes/fieldtrip-20161009/');
ft_defaults;

%%
clc;
clear all;
close all;

%%
if isempty(gcp('nocreate'))
    parpool(4,'SpmdEnabled',false);
end;

%%
rpath = '/media/samba_share/RDS/iEEG_DATA/MICRO/P05/fvSpEM/2016-11-27_11-19-21/';
files1 = dir([rpath,'lfp_dat/*.mat']);
files2 = dir([rpath,'log_dat/*Log*.mat']);

%%
load([rpath,'log_dat/',files2(2).name]);

dum = [];
dum(1,:) = LogDat1.idx(1,:);
for it = 2:size(LogDat1.idx,1)
    dum(it,:) = [LogDat1.idx(it,:)+max(LogDat1.idx(it-1,:))];
end;

if diff(LogDat1.idx,[],2) ~= diff(dum,[],2)
    error('trial assignment must match');
end;

dum2 = cell(1,size(dum,1));
for it = 1:size(dum,1)
    dum2{it} = dum(it,1):dum(it,2);
end;

tID = [dum2{:}];

tID_r = tID(ix{4});
tID_f = sort([tID(ix{5}) tID(ix{6})]);

%%
lfp_dat = cell(1,length(files1)-1);

parfor gt = 2:length(files1)
    
    fprintf([num2str(gt),'/',num2str(length(files1))]);
    dat = load([rpath,'lfp_dat/',files1(gt).name]);
    
    lfp_dat{gt} = dat.save_data{1}{1}{1};   
    
    fprintf('\n');
    
end;
lfp_dat(1) =[];

%%

[lfp_dat] = ft_appenddata([],lfp_dat{:});

%%
f_res = 0.01;
pf = 1/f_res;
    
cfg             = [];
cfg.latency     = [-1 -1e-3];
cfg.trials      = sort([tID_r]);

[dum1] = ft_selectdata(cfg,lfp_dat);

cfg             = [];
cfg.latency     = [2+1e-3 3];
cfg.trials      = sort([tID_r]);

[dum2] = ft_selectdata(cfg,lfp_dat);

cfg                     = [];
cfg.detrend             = 'yes';
cfg.demean              = 'yes';
%cfg.derivative          = 'yes';

[dum1] = ft_preprocessing( cfg , dum1 );

cfg                     = [];
cfg.detrend             = 'yes';
cfg.demean              = 'yes';
%cfg.derivative          = 'yes';

[dum2] = ft_preprocessing( cfg , dum2 );

dt = length(dum2.trial{1}) - length(dum1.trial{1});

if dt ~=0
    zp = zeros(length(dum1.label),dt/2);
    
    for it = 1:length( dum1.trial )
        dum1.trial{it} = [zp dum1.trial{it} zp];
        dum1.time{it} = 0:1/dum1.fsample:(length(dum1.trial{it})-1)/dum1.fsample;
    end;
end;

zp = zeros(length(dum1.label),length(dum1.time{1})+length(dum1.time{1})*ceil(pf/2));

for it = 1:length( dum1.trial )
    dum1.trial{it} = [zp dum1.trial{it} zp];
    dum1.time{it} = 0:1/dum1.fsample:(length(dum1.trial{it})-1)/dum1.fsample;
end;

zp = zeros(length(dum2.label),length(dum2.time{1})+length(dum2.time{1})*ceil(pf/2));

for it = 1:length( dum2.trial )
    dum2.trial{it} = [zp dum2.trial{it} zp];
    dum2.time{it} = 0:1/dum2.fsample:(length(dum2.trial{it})-1)/dum2.fsample;
end;

w = length(dum1.time{1})/dum1.fsample;

cfg                     = [];
cfg.method              = 'mtmfft';
cfg.pad                 = 'maxperlen';
cfg.taper               = 'dpss';
cfg.tapsmofrq           = 1/w;
cfg.output              = 'fourier';

[freq_b]   = ft_freqanalysis( cfg , dum1 );
[freq_ps]  = ft_freqanalysis( cfg , dum2 );

clear dum*;

cfg                     = [];
cfg.frequency           = [1 20];

[freq_b]  = ft_selectdata( cfg ,  freq_b );
[freq_ps] = ft_selectdata( cfg ,  freq_ps );

cfg                     = [];
cfg.keeptrials          = 'yes';

[pow_b]  = ft_freqdescriptives( cfg , freq_b );
[pow_ps] = ft_freqdescriptives( cfg , freq_ps );

clear freq_*;
    
%%
cfg                 = [];
cfg.variance        = 'yes';
cfg.jacknife        = 'yes';
cfg.keeptrials      = 'no';

[pow_b2] = ft_freqdescriptives( cfg , pow_b );
[pow_ps2] = ft_freqdescriptives( cfg , pow_ps );

%%
n1 = size(pow_ps.powspctrm,1);
n2 = size(pow_b.powspctrm,1);

cfg                     = [];
cfg.method              = 'montecarlo';
cfg.statistic           = 'ft_statfun_depsamplesT';
cfg.numrandomization    = 2000;
cfg.correctm            = 'none';
cfg.alpha               = 1e-3;
cfg.tail                = 0;
cfg.ivar                = 1;
cfg.uvar                = 2;

cfg.design = [];
cfg.design(1,:) = [ones(1,n1) 2*ones(1,n2)];
cfg.design(2,:) = [1:n1 1:n2];

[freq_stat] = ft_freqstatistics( cfg , pow_ps, pow_b );


%%
figure;
lm = [];
for it = 1:length(freq_stat.label)
    
    subplot(6,8,it)    
    a(it) = gca;    
    hold on;
    Y = freq_stat.mask(it,:);
    freqIdx = find(Y~=0);
    
    if sum(Y) ~=0                
        
        plot(pow_ps2.freq,pow_ps2.powspctrm(it,:),'c');
        plot(pow_ps2.freq,pow_ps2.powspctrm(it,:)+pow_ps2.powspctrmsem(it,:),'b');
        plot(pow_ps2.freq,pow_ps2.powspctrm(it,:)-pow_ps2.powspctrmsem(it,:),'b');
        
        plot(pow_ps2.freq,pow_b2.powspctrm(it,:),'g');
        plot(pow_ps2.freq,pow_b2.powspctrm(it,:)+pow_b2.powspctrmsem(it,:),'Color',[0 .25 0]);
        plot(pow_ps2.freq,pow_b2.powspctrm(it,:)-pow_b2.powspctrmsem(it,:),'Color',[0 .25 0]);
    end;
end;
axis(a,'tight');
%set(a,'YLim',[min(min(lm)) max(max(lm))]);
