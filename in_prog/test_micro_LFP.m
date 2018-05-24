%% set the path environment
addpath('/home/rouxf/tbx/fieldtrip-20170618/');
ft_defaults;

addpath('/home/rouxf/prj/Bham/code/mcode/project_EM/');
addpath('/home/rouxf/prj/Bham/code/mcode/utils/LFP/');
addpath('/home/rouxf/prj/Bham/code/mcode/visualize/continuous/');

%%
p2micDat = '/media/rouxf/rds-share/iEEG_DATA/MICRO/P07/cnEM/2017-05-04_18-27-03/lfp_dat/';
micDat = 'P07_cnEM_CSC_postHippR_2017-05-04_18-27-03_preprocANDrejectedIED_stimlocked.mat';
load([p2micDat,micDat]);

micFiles = dir([p2micDat,'*_lfp_data_',LFPdatMW.label{:},'*_stimlocked.mat']);

p2logDat = '/media/rouxf/rds-share/iEEG_DATA/MICRO/P07/cnEM/2017-05-04_18-27-03/log_dat/';
logDat = 'P07_cn_05_2017-May-04_18_24_5281500_LogFile_EMtask_LogDat.mat';
load([p2logDat,logDat]);

%%
MWdat = cell(1,length(micFiles));
for it = 1:length( micFiles )
    dat = load([p2micDat,micFiles(it).name]);
    MWdat{it} = dat.save_data{1}{1}{1};
end;

[MWdat] = ft_appenddata([],MWdat{:});

%%
[trlENC] = extract_Enc_trl(LogDat1);

dum = intersect(trlENC,selIx);
selIx = find(ismember(selIx,dum));

dum = intersect(trlENC,delIx);
delIx = find(ismember(delIx,dum));

%%
cfg                     = [];
cfg.trials              = trlENC;

[MWdat] = ft_selectdata( cfg , MWdat );

cfg                      = [];
cfg.trials               = setdiff( 1:length(MWdat.trial), delIx );

[MWdat] = ft_selectdata( cfg , MWdat );

%%
cfg                     = [];
cfg.hilbert             = 'abs';

[RMS] = ft_preprocessing( cfg, MWdat );

cfg                     = [];
cfg.keeptrials          = 'no';

[RMS] = ft_timelockanalysis( cfg , RMS );

cfg                     = [];
cfg.keeptrials          = 'no';
cfg.covariance          = 'yes';

[COV] = ft_timelockanalysis( cfg , MWdat );

%%
BipRefMWdat = MWdat;
BipRefMWdat.trial = {};
BipRefMWdat.label = {};
NchanPairs = length( MWdat.label )-1;
Nsamp = length( MWdat.time{1} );

for jt = 1:length( MWdat.trial)
    
    fprintf([num2str(jt),'/',num2str(length( MWdat.trial))]);
    x = MWdat.trial{jt};
    
    BipRefMWdat.trial{jt} = zeros(NchanPairs,Nsamp);
    
    for it = 1:NchanPairs
        
        s1 = x(it,:);
        s2 = x(it+1,:);
        
        BipRefMWdat.trial{jt}(it,:) = s1-s2;
        clear s1 s2;
        
    end;
    clear x;
    
    
    fprintf('\n');
    
end;
for it = 1:NchanPairs
    BipRefMWdat.label(it) = {[MWdat.label{it},'-',MWdat.label{it+1}]};
end;

%%
cfg                     =[];
cfg.trials              = selIx;

[LFPdatMW] = ft_selectdata( cfg , LFPdatMW );

%%
[erpDat1,TFR1] = compute_ERP_and_TFR_macroData_EM(LFPdatMW);
[erpDat2,TFR2] = compute_ERP_and_TFR_macroData_EM(BipRefMWdat);

cfg                     = [];
cfg.output              = 'pow';
cfg.keeptrials          = 'yes';
cfg.toi                 = [-1 4.5];

[TFR1,~]=sh_subtr1of(cfg,TFR1);

[TFR2,~]=sh_subtr1of(cfg,TFR2);

%%
cfg                     = [];
cfg.latency             = [2 max(LFPdatMW.time{1})];
[dum] = ft_selectdata( cfg, LFPdatMW );

dt = 250;

cfg                     = [];
cfg.bpfilter            ='yes';
cfg.bpfreq              = [50 70];

[filt] = ft_preprocessing( cfg , dum );

epch2 = [];
c = 0;
for it = 1:length( filt.trial )
    fprintf([num2str(it),'/',num2str(length(filt.trial))]);
    
    [mix1] = local_max(filt.trial{it});
    [mix2] = local_max(-filt.trial{it});
    
    for jt = 1:length(mix1)
        selIx = mix1(it)-dt:mix1(it)+dt;
        
        if (min(selIx)>0) && (max(selIx) < length(dum.trial{it}))
            c = c + 1;
            epch2(c,:) = dum.trial{it}(selIx);
        end;
        
    end;
    
    fprintf('\n');
end;

%%
dum2 = struct;
dum2.label = {'sta_chan_post'};
dum2.fsample = 1e3;
dum2.cfg = [];
for it = 1:size( epch2 ,1)
    dum2.trial{it} = [epch2(it,:) epch2(it,:) epch2(it,:) epch2(it,:) epch2(it,:)];
    dum2.time{it} = 0:1/dum2.fsample:(length(dum2.trial{it})-1)/dum2.fsample;
end;

%%
[erpDat2,TFR2] = compute_ERP_and_TFR_macroData_EM(dum2);

%%
cfg                     = [];
cfg.output              = 'pow';
cfg.keeptrials          = 'yes';
cfg.toi                 = [min(TFR2.time) max(TFR2.time)];

[TFR2,slp]=sh_subtr1of(cfg,TFR2);

%%
selIx = size(epch2,2)*1+1:size(epch2,2)*2;

cfg                     = [];
cfg.latency             = [min(dum2.time{1}(selIx)) max(dum2.time{1}(selIx))];

[TFR2] = ft_selectdata( cfg , TFR2 );

%%