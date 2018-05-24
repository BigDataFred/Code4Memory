%%

addpath('/media/rouxf/rds-share/Common/fieldtrip-20170115/');
ft_defaults;

%% subject ID
pID = {'P02'};

%% path to main subjet folder
rpath = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',pID{1},filesep];

[exp,sesh] = get_EM_rec_info(rpath);%extract subfolders from main folder

[rpath] = [rpath,exp,filesep]; % set root path for data reading

gt = 2;% select a session

%% load the LOG-data
path_2_log = [rpath,sesh{gt},filesep,'log_dat',filesep];
fn = dir([path_2_log,'*_EMtask_LogDat.mat']);
[LogDat] = load([path_2_log,fn.name]);%load the data

%% do some processing on the LOG-data
for it = 1:length(LogDat.ix)
    if size(LogDat.ix{it},2) > size(LogDat.ix{it},1)
        LogDat.ix{it} = LogDat.ix{it}';
    end;
end;

trlIx1 = [LogDat.ix{1};LogDat.ix{2};LogDat.ix{3}];% trial indexes for cond1-3
trlIx2 = [LogDat.ix{4};LogDat.ix{5};LogDat.ix{6}];% trial indexes for hits & misses

% extract the trial information for the encoding phase
[trl_all] = extract_Enc_trl(LogDat.LogDat1);% all, remembered (r), not-remembered (nr)

%% SPK-data source path and data files
path_2_spk = [rpath,sesh{gt},filesep,'spike_dat',filesep];
spk_files = dir([path_2_spk,'*.mat']);

%% LFP-data source path and data files
path_2_lfp = [rpath,sesh{gt},filesep,'lfp_dat',filesep];
lfp_files = dir([path_2_lfp,'*.mat']);

%% load the CLUSTER-INFO
p2d = rpath;
load([p2d,'Cluster_info_EM_',pID{1},'_',sesh{gt},'.mat']);

%% load the LFP-data
lfp_dat = cell(1,length(lfp_files));
for ft = 1:length(lfp_files  )
    fprintf([num2str(ft),'/',num2str(length(lfp_files))]);
    load([path_2_lfp,lfp_files(ft).name]);
    lfp_dat{ft} = save_data{1}{1}{1};
    clear save_dat;
    fprintf('\n');
end;

[lfp_dat] = ft_appenddata([],lfp_dat{:});

%% do some processing on the LFP-data
cfg                     = [];
cfg.trials              = trl_all;

[lfp_dat] = ft_selectdata( cfg, lfp_dat );

% cfg                     = [];
% cfg.bpfilter            = 'yes';
% cfg.bpfreq              = [1 100];
% 
% [lfp_dat] = ft_preprocessing( cfg , lfp_dat );

%%
% cfg                     = [];
% 
% [cfg] = ft_databrowser(cfg,lfp_dat);
% 
% [lfp_dat] = ft_rejectartifact(cfg, lfp_dat );
%%
cfg                     = [];
cfg.keeptrials          = 'yes';

[tlck] = ft_timelockanalysis( cfg , lfp_dat );

%%
cfg                     = [];
cfg.method              = 'mtmconvol';
cfg.output              = 'pow';
cfg.pad                 = 'maxperlen';
cfg.taper               = 'dpss';
cfg.tapsmofrq           = 1;
cfg.foi                 = [1:20];
cfg.toi                 = lfp_dat.time{1}(1):0.067:lfp_dat.time{1}(end);
cfg.t_ftimwin           = ones(1,length(cfg.foi));
cfg.keeptrials          = 'yes';

[TFR1] = ft_freqanalysis( cfg , lfp_dat );

cfg                     = []; 
cfg.channel             = {'all','-SYNC'};

[TFR1] = ft_selectdata( cfg , TFR1 );

cfg                     = [];
cfg.baseline            = [-1 0];
cfg.baselinetype        = 'relchange';

[bTFR1] = ft_freqbaseline( cfg , TFR1 );

%%
cfg                     = [];
cfg.method              = 'mtmconvol';
cfg.output              = 'pow';
cfg.pad                 = 'maxperlen';
cfg.taper               = 'dpss';
cfg.tapsmofrq           = 10;
cfg.foi                 = [20:4:100];
cfg.toi                 = lfp_dat.time{1}(1):0.067:lfp_dat.time{1}(end);
cfg.t_ftimwin           = 0.25*ones(1,length(cfg.foi));
cfg.keeptrials          = 'yes';

[TFR2] = ft_freqanalysis( cfg , lfp_dat );

cfg                     = []; 
cfg.channel             = {'all','-SYNC'};

[TFR2] = ft_selectdata( cfg , TFR2);

cfg                     = [];
cfg.baseline            = [-1 0];
cfg.baselinetype        = 'relchange';

[bTFR2] = ft_freqbaseline( cfg , TFR2 );

%% load the SPK-data
spk_dat = cell(1,length(spk_files));
for ft = 1:length(spk_files  )
    fprintf([num2str(ft),'/',num2str(length(spk_files))]);
    load([path_2_spk,spk_files(ft).name]);
    spk_dat{ft} = save_data{1}{1}{1};   
    clear save_dat;
    fprintf('\n');
end;

%% indexes to group microwires into benke-fried banks
selIx = 1:8:length(spk_files);
bfIdx = {};
for it = 1:length( selIx )
    bfIdx{it} = [selIx(it):selIx(it)+7];% always 8 MWs belong to 1 BF bank
end;

%% select the SU 

SUidx = intersect(find(clusterInfo(:,5)==1),find(clusterInfo(:,9)==1));% SUs must have SU properties and also show increase in response
%SUidx = find(clusterInfo(:,9)==1);% SUs must have SU properties and also show increase in response
%SUidx = SUidx(6);
selIx = [clusterInfo(SUidx,2) clusterInfo(SUidx,3)+1];% indexes for each cluster per MW

%% compute the spike triggered average (STA)
itc = {};
STA = {};

mITC = [];
for ct = [3 4]%1:size(selIx,1)
    
    %%
    bfSELIdx = [];
    for it = 1:length(bfIdx)
        bfSELIdx(it) = ismember(selIx(ct,1),bfIdx{it});
    end;
    bfSELIdx = find(bfSELIdx);
    
    %%
    tS = spk_dat{selIx(ct,1)}.time{selIx(ct,2)}.*1e3;
    trl = spk_dat{selIx(ct,1)}.trial{selIx(ct,2)};
    trlID = 1:length(trl_all);
    
    spkRaster = [];
    dt = -2000:5000;
    for kt = 1:length( trlID )
        
        ix = find(trl == trlID(kt));
        x = tS(ix);
        x(x<dt(1)) = [];
        x(x>dt(end)) = [];
        
        spkRaster(kt,:) = hist(x,dt);
        
    end;
    
    %%
    dt2 = 1/4*3*lfp_dat.fsample;
    c=0;
    nTRL = size(spkRaster,1);
    for kt = 1:size(spkRaster,1)
        fprintf([num2str(kt),'/',num2str(size(spkRaster,1))]);
        ix = [];
        ix = find(spkRaster(kt,:));
        
        if ~isempty(ix)
            ix = ix(dt(ix)>2);
            for nt = 1:length(ix)
                start = ix(nt)-dt2;
                stop = ix(nt)+dt2;
                if (start >1) && (stop < size(lfp_dat.trial{kt},2))
                    c = c+1;
                    STA{ct}(c,:,:) = lfp_dat.trial{kt}(bfIdx{bfSELIdx},start:stop);
                end;
            end;
        end;
        fprintf('\n');
    end;
        
    %%
    h = zeros(size(STA{ct}));
    for it = 1:size(STA{ct},1)
        for jt = 1:size(STA{ct},2)
            h(it,jt,:) = hilbert(squeeze(STA{ct}(it,jt,:)));           
        end;
    end;
    h = h./abs(h);
    itc{ct} = squeeze(abs(sum(h,1))./size(h,1));   
    
end;

%%
figure;
for it = 1:length(bfIdx)-1
    subplot(length(bfIdx)-1,1,it);
    imagesc(TFR.time,TFR.freq,squeeze(median((bTFR.powspctrm(bfIdx{it},:,:)),1)));
    axis xy;
end;

%%
ct = 3;
dum = mean(itc{ct},2);
[i1] = find(dum == max(dum));
mITC(ct) = dum(i1);
ix = find(selIx(ct,1)==bfIdx{bfSELIdx});
ix = 2;

figure;
subplot(121);
hold on;
plot(-dt2:dt2,squeeze(mean(STA{ct},1)),'k');
plot(-dt2:dt2,mean(squeeze(STA{ct}(:,ix,:)),1),'r');
axis tight;
subplot(122);
plot(-dt2:dt2,mean(squeeze(STA{ct}(:,i1,:)),1),'r');
axis tight;

figure;
imagesc(-dt2:dt2,1:size(STA{ct},1),squeeze(STA{ct}(:,ix,:)));
axis xy;