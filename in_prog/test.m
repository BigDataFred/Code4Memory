%% set the path def
restoredefaultpath;
addpath(genpath('C:\toolbox\AnalysisFred\'));% custom code
addpath(genpath('C:\toolbox\Nlx\MatlabImportExport_v6.0.0\'));%needed to read Nlx data
addpath(genpath('C:\toolbox\osort-v3.0\osort-v3-rel\code\'));%needed for spike detection & sorting
addpath(genpath('C:\toolbox\chronux_2_12\'));% needed for spectral analysis
addpath(genpath('C:\toolbox\wave_clus\Wave_clus\'));%needed for spike detection & sorting
addpath('C:\toolbox\fieldtrip-20160309\');% spectral analysis, signal processing, spike detection
ft_defaults;

%% open pool of workers
%parpool('SpmdEnabled',false);
clear;
clc;

%% set the data-paths directions
p2Nlxdata = 'C:\Data\Nlx\tuning\';
p2logdat = 'C:\Experiments\EMpairs_v4_2016-1007\log\Tuning\';
p2logparams = 'C:\Experiments\EMpairs_v4_2016-1007\params\';

%% select session to analyze
sesh = 1;
mode = 'stimlocked';%'resplocked';
pre = 1.5;
post = 2.5;
chanSel = [ 1:48 ];

makeLFP = 1;

%% create the session labels of the Nlx data
Nlxdat = {};
Nlxdat{1} = '2016-07-09_11-13-28\';
Nlxdat{2} = 'P02_TS01_2016-Jul-10_10-50-00\';
Nlxdat{3} = 'P02_TS07_2016-Jul-18_11_56_08\';
Nlxdat{4} = 'P02_TS08_2016-Jul-19_10_54_39\';

%% create the session labels of the log data
Tlogfile = {};
Tlogfile{1} = 'P02_TS01_log_ctune_09-Jul-2016_11_40_xx.txt';
Tlogfile{2} = 'P02_TS02_log_ctune_10-Jul-2016_10_46_xx.txt';
Tlogfile{3} = 'P02_TS07__log_ctune_18-Jul-2016_11_56_08.txt';
Tlogfile{4} = 'P02_TS09__log_ctune_19-Jul-2016_10_54_39.txt';

%% get the logfile data
f.p2logf = 'C:\Experiments\EMpairs_v4_2016-1007\log\Tuning\';
f.logf = Tlogfile{sesh};
p.ncols = 8;

[LogDat] = getNewLogDataTune(f,p);

if length(unique(LogDat.dat(find(LogDat.stimID == LogDat.ID(1)),2))) >1
    error('wrong stimulus event assignment');
end;

%% get the event file data
[event] = ft_read_event([p2Nlxdata,Nlxdat{sesh}]);

ttl_idx = find([event(:).value] == 7);% get onset events

if LogDat.n ~= length(ttl_idx)
    error('number of events is out of range');
end;

%% find the offset events
k = 0;
ttl_idx2 = zeros(length(ttl_idx),1);
for it = 1:length(ttl_idx)
    for jt = 1:2
        if event(ttl_idx(it)+jt).value ==0
            k = k+1;
            ttl_idx2(k) = ttl_idx(it)+jt;
            break;
        end;
    end;
end;

if length(ttl_idx2) ~= length(ttl_idx)
    error('number of events must be equal');
end;

%% extract Fs and ADV
[hdr] = ft_read_header([p2Nlxdata,Nlxdat{sesh}]);

%% convert event ts to samples
for it = 1:length([event(:).value])
    event(it).samples = (event(it).timestamp - hdr.FirstTimeStamp)/hdr.TimeStampPerSample + 1;
end;

if unique([event(ttl_idx).value]) ~=7
    error('trigger values are out of range');
end;

if length([event(ttl_idx).value]) ~=LogDat.n
    error('number of events out of range');
end;

%%
switch mode
    case 'stimlocked'
        trl = zeros(length(ttl_idx),3);
        for it =1:length(ttl_idx)
            trl(it,:) = [event(ttl_idx(it)).sample-(2*pre*hdr.Fs) event(ttl_idx(it)).sample+(2*pre*hdr.Fs) -(2*pre*hdr.Fs)];
        end;
        del_idx = [find(sign(trl(:,1))==-1) find(trl(:,2)>event(end).sample+2*post*hdr.Fs)];
        trl(del_idx,:) = [];
        LogDat.RT(del_idx) = [];
        %%
    case 'resplocked'
        trl = zeros(length(ttl_idx2),3);
        for it =1:length(ttl_idx2)
            trl(it,:) = [event(ttl_idx2(it)).sample-(2*pre*hdr.Fs) event(ttl_idx2(it)).sample+(2*pre*hdr.Fs) -(2*pre*hdr.Fs)];
        end;
        del_idx = [find(sign(trl(:,1))==-1) find(trl(:,2)>event(end).sample+2*post*hdr.Fs)];
        trl(del_idx,:) = [];
        LogDat.RT(del_idx) = [];
        ttl_idx = ttl_idx2;
end;

%% filter settings
Wp = [600 8000]./hdr.Fs/2;
[b,a] = butter(4,Wp);

Hd{1} = b;
Hd{2} = a;

%%
files = dir([p2Nlxdata,Nlxdat{sesh},'*.ncs']);
%parpool;
tlck = cell(1,length(chanSel));
%%
for ct = 1:length(chanSel)
    
    files;
    Nlxdat;
    LogDat2 = LogDat;
    trl2 = trl;
    ttl_idxc = ttl_idx;
    
    %% st = tic;
    fprintf([num2str(ct),'/',num2str(length(chanSel))]);
    fprintf('\n');
    
    csc_filename = files(chanSel(ct)).name;
    
    
    if makeLFP ==1
        %% make the LFP data
        cfg             = [];
        cfg.dataset     = [p2Nlxdata,Nlxdat{sesh},csc_filename];
        %cfg.dftfilter   = 'yes';
        %cfg.dftfreq     = [50:50:300];
        cfg.padtype     = 'data';
        cfg.padding     = 100;
        
        [dum]           = ft_preprocessing(cfg);
        
        %% high-pass
        cfg = [];
        cfg.hpfilter        = 'yes';
        cfg.hpfreq          = 1/30;
        cfg.hpfiltord       = 2;
        cfg.padtype         = 'mirror';
        cfg.padding         = 100;
        
        [dum]               = ft_preprocessing(cfg,dum);
        
        %% low-pass
        cfg = [];
        cfg.lpfilter        = 'yes';
        cfg.lpfreq          = 200;
        cfg.lpfiltord       = 2;
        cfg.padtype         = 'mirror';
        cfg.padding         = 100;
        
        [dum]               = ft_preprocessing(cfg,dum);
        
        %% trial segmentation
        cfg             = [];
        cfg.trl         = trl2;
        
        [lfp]           = ft_redefinetrial(cfg,dum);
        %clear dum;
        
        %% downsample LFP
        cfg             = [];
        cfg.resamplefs  = 1200;
        
        [lfp]           = ft_resampledata(cfg,lfp);
        
        %%
        cfg             = [];
        cfg.keeptrials  = 'yes';
        
        [dum]          = ft_timelockanalysis(cfg,lfp);
        
        %% find outlier trials
        x = squeeze(mean(abs(dum.trial).^2,3));
        z = (x - mean(x))./std(x);
        del_idx = find(abs(z) >= 2.5);
        del_idx = [];
        %% remove outlier trials from data
        trl2(del_idx,:) = [];
        
        ttl_idxc(del_idx) = [];
        
        LogDat2.RT(del_idx) = [];
        LogDat2.stimID(del_idx) = [];
        LogDat2.dat(del_idx,:) = [];
        
        cfg             = [];
        cfg.trials      = setdiff(1:size(dum.trial,1),del_idx);
        
        [lfp]           = ft_redefinetrial(cfg,lfp);
        
        %% sort data according to RTs
        [~,s_idx] = sort(LogDat2.RT);
        
        trl2 = trl2(s_idx,:);
        
        ttl_idxc = ttl_idxc(s_idx);
        
        lfp.trial = lfp.trial(s_idx);
        
        LogDat2.RT = LogDat2.RT(s_idx);
        LogDat2.stimID = LogDat2.stimID(s_idx);
        LogDat2.dat = LogDat2.dat(s_idx,:);
        
        if length(unique(LogDat2.dat(find(LogDat2.stimID == LogDat2.ID(1)),2))) >1
            error('wrong stimulus event assignment');
        end;
        
        %% make the evoked potential
        cfg             = [];
        cfg.keeptrials  = 'yes';
        
        [tlck{ct}]          = ft_timelockanalysis(cfg,lfp);
        
        visualize_RTvsLFP(tlck{ct},LogDat2);
    end;
end;