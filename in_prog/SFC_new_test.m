restoredefaultpath;

fn = dir('/media/rouxf/rds-share/Common/fieldtrip-*');
addpath(['/media/rouxf/rds-share/Common/',fn.name]);% spectral analysis, signal processing, spike detection
ft_defaults;
addpath(genpath('/media/rouxf/rds-share/Fred/code/mcode/custom/'));

%%
xt = 1;
pID = {'P23AMS'};

rpath = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',pID{1},filesep];

[exp,sesh] = get_EM_rec_info(rpath);

[rpath] = [rpath,exp,filesep];

p2LFP = [rpath,sesh{xt},'/lfp_dat/'];
LFPfiles = dir([p2LFP,'*.mat']);

p2SPK = [rpath,sesh{xt},'/spike_dat/'];
SPKfiles = dir([p2SPK,'*.mat']);


[savepath] = rpath;

%% load the Log data
p2LOG = [rpath,sesh{xt},filesep,'log_dat',filesep];
fn = dir([p2LOG,'*_EMtask_LogDat.mat']);
[LogDat] = load([p2LOG,fn.name]);%load the data

%% do some processing on the Log data
for it = 1:length(LogDat.ix)
    if size(LogDat.ix{it},2) > size(LogDat.ix{it},1)
        LogDat.ix{it} = LogDat.ix{it}';
    end;
end;

% extract the trial information for the encoding phase
[trlAll] = extract_Enc_trl(LogDat.LogDat1);% all, remembered (r), not-remembered (nr)

%% load the LFP data
LFPdat = cell(1,length(LFPfiles));
parfor it = 1:length(LFPfiles)
    fprintf([num2str(it),'/',num2str(length(LFPfiles))]);
    dat = load([p2LFP,LFPfiles(it).name]);
    LFPdat{it} = dat.save_data{1}{1}{1};
    fprintf('\n');
end;

[LFPdat] = ft_appenddata( [] , LFPdat{:} );

%% append the LFP channels
cfg                     = [];
cfg.trials              = trlAll;

[LFPdat] = ft_selectdata( cfg , LFPdat );

%% load the SPK data
SPKdat = cell(1,length(SPKfiles));
parfor it = 1:length(SPKfiles)
    fprintf([num2str(it),'/',num2str(length(SPKfiles))]);
    dat = load([p2SPK,SPKfiles(it).name]);
    SPKdat{it} = dat.save_data{1}{1}{1};
    fprintf('\n');
end;

%%
for it = 1:length( SPKdat )
    if ~isempty([SPKdat{it}.timestamp{:}])
        cfg                     = [];
        cfg.trials              = trlAll;
        
        [SPKdat{it}] = ft_spike_select( cfg , SPKdat{it} );
    end;
end;

%%
load([rpath,'Cluster_info_EM_',pID{1},'_',sesh{xt},'.mat']);

SUidx = intersect(find(clusterInfo(:,5)==3),find(clusterInfo(:,14)==1));

%ChanIx = find(strcmp(LFPdat.label,'CSC_RA3'));
%SUidx = find(clusterInfo(:,2)==ChanIx);

selIx = [clusterInfo(SUidx,2) clusterInfo(SUidx,3)+1];
%selIx = selIx([3 6],:);

%%
statStsL = cell(1,size( selIx ,1 ));
statStsH = cell(1,size( selIx ,1 ));
staPre  = cell(1,size( selIx ,1 ));
staPost = cell(1,size( selIx ,1 ));
CSCix   = cell(1,size( selIx ,1 ));

STA = cell(1,size( selIx ,1 ));
itc   = cell(1,size( selIx ,1 ));

%%
for ut = 1:size( selIx ,1)
    
    cfg                     = [];
    cfg.spikechannel        = SPKdat{selIx(ut,1)}.label(selIx(ut,2));
    
    dum = ft_spike_select( cfg , SPKdat{selIx(ut,1)} );
    
    %%
    source = dum.hdr.label;
    source(regexp(source,'\d{1,3}')) = '*';
    BFchan = ft_channelselection({source},LFPdat.label);
    
    cfg                     = [];
    cfg.channel             = BFchan;
    
    [dum2] = ft_selectdata( cfg , LFPdat );
    
    %%
    dt = min(LFPdat.time{1}):1/LFPdat.fsample:max(LFPdat.time{1});
    spk = cell(1,length(LFPdat.trial));
    for it = 1:length(trlAll)
        
        ix = find(dum.trial{1} == trlAll(it));
        x = dum.time{1}(ix);
        
        n = hist(x,dt);
        
        spk{it} = n;
        unique(spk{it})
    end;
    
    %%
    nsamp = length(dum2.time{1});
    for kt = 1:length( spk )
        
        ix = find(spk{kt});
        ix = ix';
        
        trl = [ix-1000 ix+1000 -1000*ones(length(ix),1)];
        del1 = find(sign(trl(:,1))==-1);
        del2 = find(trl(:,2)>nsamp);
        trl([del1;del2],:) = [];
        
        if ~isempty(trl)
                                    
%             cfg                     = [];
%             cfg.trl                 = trl;
%                         
%             [dum3] = ft_redefinetrial( cfg , dum2 );
            
            dum3 = dum2;
            dum3.trial ={};
            dum3.time  ={};
            for mt = 1:size(trl,1)
                dum3.trial{mt} = dum2.trial{kt}(:,trl(mt,1):trl(mt,2));
                dum3.time{mt} = -1:1/dum2.fsample:1;
            end;
            
            if kt >1
                [sta] = ft_appenddata([],sta,dum3);
            else
                sta = dum3;
            end;
        end;
        
    end;
    
    df = dum2.fsample/2/1024;
    
    cfg                     = [];
    cfg.method              = 'wavelet';
    cfg.output              = 'fourier';
    cfg.keeptrials          ='yes';
    cfg.foi                 = df:df:45;
    cfg.toi                 =-1:0.01:1;    
    cfg.width               = 5;
    cfg.pad                 = 'nextpow2';
    
    [freq] = ft_freqanalysis(cfg , sta );
    
    itc{ut}                 = [];
    itc{ut}.freq            = freq.freq;
    itc{ut}.label           = freq.label;
    itc{ut}.time            = freq.time;
    itc{ut}.dimord          = 'chan_freq_time';
    itc{ut}.itc             = freq.fourierspctrm;
    itc{ut}.itc             = itc{ut}.itc./abs(itc{ut}.itc);
    itc{ut}.itc             = abs(mean(itc{ut}.itc,1));
        
    cfg                     = [];
    cfg.keeptrials          = 'yes';
    
    [STA{ut}] = ft_timelockanalysis( cfg , sta );
    
%     %%
%     LFPdatAll = dum2;
%     LFPdatAll.label(end+1) = dum.label;
%     
%     for it = 1:length( spk )
%         
%         LFPdatAll.trial{it} = [LFPdatAll.trial{it};spk{it}];
%         
%     end;
%     
% %     dum2.hdr = dum.hdr;
% %     trl = dum2.cfg.previous.previous.previous{1}.previous.previous.trl;
% %     trl = trl(trlAll,:);
% %     dum2.cfg.previous.previous.previous{1}.previous.previous.trl = trl;
% %     
% %     [LFPdatSPK] = ft_appendspike([]  ,dum2 , dum);
%     
%     %%
%     cfg                     = [];
%     cfg.timwin              = [-0.25 0.25];
%     cfg.spikechannel        = dum.label;
%     cfg.channel             = LFPdatAll.label(1:end-1);
%     cfg.latency             = [-1 0];
%     cfg.keeptrials         = 'yes';
%     
%     [staPre{ut}] = ft_spiketriggeredaverage( cfg , LFPdatAll );
%     
%     cfg                     = [];
%     cfg.timwin              = [-0.25 0.25];
%     cfg.spikechannel        = dum.label;
%     cfg.channel             = LFPdatAll.label(1:end-1);
%     cfg.latency             = [2 5];
%     cfg.keeptrials         = 'yes';
%     
%     [staPost{ut}] = ft_spiketriggeredaverage( cfg , LFPdatAll );
%     
%     CSCix{ut} = find(strcmp(staPost{ut}.label,dum.hdr.label));
%     
%      %%
%     cfg                     = [];
%     cfg.spikechannel        = LFPdat.label;
%     cfg.method              = 'linear';
%     cfg.timwin              = [-0.0005 0.0005];
%     cfg.spikechannel        = dum.label;
%     cfg.channel             = LFPdatAll(1:end-1);
%     
%     [LFPdatAll] = ft_spiketriggeredinterpolation( cfg , LFPdatAll );
%     
%     for it = 1:length( spk )
%         
%         LFPdatAll.trial{it}(end,:) = spk{it};
%         
%     end;
%     
% %     %%
% %     for it = 1:length(LFPdatAll.trial)
% %         chck = unique(LFPdatAll.trial{it}(end,:));
% %         if any(chck >0) && any(chck<1)
% %             ix = find(LFPdatAll.trial{it}(end,:) > 0 & LFPdatAll.trial{it}(end,:) < 1);
% %             LFPdatAll.trial{it}(end,ix) = 1;
% %         end;
% %     end;
% %     
% %     for it = 1:length(LFPdatAll.trial)
% %         chck = unique(LFPdatAll.trial{it}(end,:));
% %         
% %     end;
%     
%     %%
%     cfg                     = [];
%     cfg.latency             = [-1 0];
%     
%     [preSTIM] = ft_selectdata( cfg , LFPdatAll );
%     
%     cfg                     = [];
%     cfg.latency             = [2 5];
%     
%     [postSTIM] = ft_selectdata( cfg , LFPdatAll );
%     
%     cfg                                 = [];
%     cfg.method                          = 'mtmfft';
%     cfg.foilim                          = [2 20];
%     cfg.timwin                          = [-0.25 0.25];
%     cfg.spikechannel                    = LFPdatAll.label{end};
%     cfg.channel                         = LFPdatAll.label(1:end-1); 
%     cfg.taper                           = 'hanning';
%     
%     [stsConvol1] = ft_spiketriggeredspectrum( cfg , preSTIM );
%     [stsConvol2] = ft_spiketriggeredspectrum( cfg , postSTIM );
%     
%     %%
%     cfg                     = [];
%     cfg.spikechannel        = dum.label;
%     cfg.channel             = LFPdatAll.label(1:end-1); 
%     cfg.method              = 'plv';
%     cfg.avgoverchan         ='weighted';
%     %cfg.winstepsize         = 0.067;
%     %cfg.timwin              = 0.25;
%     
%     [statStsL{ut,1}] = ft_spiketriggeredspectrum_stat( cfg , stsConvol1 );
%     [statStsL{ut,2}] = ft_spiketriggeredspectrum_stat( cfg , stsConvol2 );      
%     
%     %%
%     cfg                     = [];
%     cfg.latency             = [-1 0];
%     
%     [preSTIM] = ft_selectdata( cfg , LFPdatAll );
%     
%     cfg                     = [];
%     cfg.latency             = [2 5];
%     
%     [postSTIM] = ft_selectdata( cfg , LFPdatAll );
%     
%     cfg                                 = [];
%     cfg.method                          = 'mtmfft';
%     cfg.foilim                          = [20 100];
%     cfg.timwin                          = [-0.05 0.05];
%     cfg.spikechannel                    = LFPdatAll.label{end};
%     cfg.channel                         = LFPdatAll.label(1:end-1); 
%     cfg.taper                           = 'hanning';
%     
%     [stsConvol1] = ft_spiketriggeredspectrum( cfg , preSTIM );
%     [stsConvol2] = ft_spiketriggeredspectrum( cfg , postSTIM );
%     
%     %%
%     cfg                     = [];
%     cfg.spikechannel        = dum.label;
%     cfg.channel             = LFPdatAll.label(1:end-1); 
%     cfg.method              = 'plv';
%     cfg.avgoverchan         ='weighted';
%     %cfg.winstepsize         = 0.067;
%     %cfg.timwin              = 0.25;
%     
%     [statStsH{ut,1}] = ft_spiketriggeredspectrum_stat( cfg , stsConvol1 );
%     [statStsH{ut,2}] = ft_spiketriggeredspectrum_stat( cfg , stsConvol2 );
    
end;

%%
for it = 1:length(STA)
    ix  = find(strcmp(SPKdat{selIx(it,1)}.hdr.label,itc{it}.label));
    for jt = 1:length(STA{it}.label)
        figure;
        subplot(4,1,1:3);
        imagesc(itc{it}.time,itc{it}.freq,squeeze(mean(itc{it}.itc(:,jt,:,:),2)));
        axis xy;
        subplot(4,1,4);
        plot(STA{it}.time,squeeze(mean(STA{it}.trial(:,jt,:),1)));
        axis tight; caxis([0 .15]);        
    end;
    figure;
    subplot(4,1,1:3);
    imagesc(STA{it}.time,1:size(STA{it}.trial,1),squeeze(STA{it}.trial(:,ix,:)));
    subplot(4,1,4);
    plot(STA{it}.time,squeeze(mean(STA{it}.trial(:,ix,:),1)));
    axis tight;
end;

%%
for it = 1:length(staPost)
    figure;
    subplot(121);
    hold on;
    plot(staPre{it}.time,staPre{it}.avg,'Color',[.75 .75 .75]);
    plot(staPre{it}.time,staPre{it}.avg(CSCix{it},:),'r');    
    axis tight;
    
    subplot(122);
    hold on;
    plot(staPost{it}.time,staPost{it}.avg,'Color',[.75 .75 .75]);
    plot(staPost{it}.time,staPost{it}.avg(CSCix{it},:),'r');    
    axis tight;
    
end;