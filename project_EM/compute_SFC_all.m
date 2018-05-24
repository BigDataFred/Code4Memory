%% reset workspace
clear all;
%close all;
clc;
date
%% set Matlab path
restoredefaultpath;
addpath(genpath('/media/rouxf/rds-share/Fred/code/mcode/custom/'));
addpath(genpath('/media/rouxf/rds-share/Fred/code/mcode/toolboxes/chronux_2_11/'));
addpath('/media/rouxf/rds-share/Fred/code/mcode/toolboxes/fieldtrip-20161009/');
ft_defaults;

%% extract the source-path information and directory names
pID = 'P04';% patient ID must be set manually
rpath = ['',filesep,'media',filesep,'rouxf',filesep,'rds-share',filesep,'iEEG_DATA',filesep,'MICRO',filesep,pID,filesep];

[exp,sesh] = get_EM_rec_info(rpath);
exp
sesh
rpath = [rpath,exp,filesep]
%%
ct = 0;%counter
CSC ={};%save labels of significant channels
Cgrm = {};

for gt = 1:length(sesh)%loop over the seshions
    
    %% spk-data source path and data files
    path_2_spk = [rpath,sesh{gt},filesep,'spike_dat',filesep];
    spk_files = dir([path_2_spk,'*.mat']);
    
    %% load the Log data
    path_2_log = [rpath,sesh{gt},filesep,'log_dat',filesep];
    fn = dir([path_2_log,'*_LogFile_EMtask_LogDat.mat']);
    
    load([path_2_log,fn.name]);%load the data
    
    %% extract the trial information for the encoding phase
    [trl_all,trl_r,trl_nr] = extract_Enc_trl(LogDat1,ix);% all, remembered (r), not-remembered (nr)
    
    %% load the LFP data for the selected seshion
    path_2_lfps = [rpath,sesh{gt},filesep,'lfp_dat',filesep];
    lfp_files = dir([path_2_lfps,'*.mat']);
    
    [lfp_dat] = load_lfp_data(lfp_files,path_2_lfps,trl_all);
    
    %%    loop over all CSC channels
    spike_dat = cell(1,length(spk_files));
    for zt = 1:length(spk_files)
        
        fprintf([num2str(zt),filesep,num2str(length(spk_files))]);
        
        % load the spike data
        [dat] = load([path_2_spk,spk_files(zt).name]);
        
        [spike_dat] = dat.save_data{1}{1}{1};% structure with spike data
        clear dat;
        spike_dat
        
        %%        
        cfg                             = [];
        cfg.channel                     = spike_dat.label;
        
        [dum] = ft_selectdata( cfg, lfp_dat );
        
        %%
        Fs = 1000; % Sampling Rate of the binned point process
        dt = -2000:1:4000;% used for the binning of spike times (PSTH)
        cID = unique(spike_dat.unit{1});% IDs of individual clusters
        
        %%       
        for it = 1:length(cID) % loop over all clusters
            
            % get the timestamps for each cluster
            i_clus = find(spike_dat.unit{1} == cID(it));
            
            % original timestamps are in micro-seconds 1e-6
            ts_all = (spike_dat.timestamp{1}( i_clus )./1e3);% save the timestamps in ms
            
            ts_trl  = cell(1,length(trl_all));
            
            c = 0;  wvf = [];
            for jt = 1:length(trl_all)%loop over the trials
                
                i_trl = find(spike_dat.trial{1}(i_clus) == trl_all(jt));% find all the timestamps for selected trial
                
                ts_trl{jt} = spike_dat.time{1}( i_clus( i_trl )).*1e3;% save the timestamps in ms
                
                ts_trl{jt}(ts_trl{jt}<dt(1)) =[];% remove timestamps that lie outside of edges
                ts_trl{jt}(ts_trl{jt}>dt(end)) =[];
                
                for kt = 1:length(i_trl)
                    c = c+1;
                    [wvf(c,:)] = squeeze(spike_dat.waveform{1}(:,:,i_clus( i_trl(kt) )))';% extract the spike waveform for each cluster
                end;
            end;
            mV = max(abs(mean(wvf,1)));
            
            %% ISI for all timestamps
            isi = diff([ts_trl{:}]);
            isi = isi(isi>0);
            isi = isi(isi<=1000);
            
            dt2 = 0:1000;%bins in ms
            [n_isi] = hist(isi,dt2);
            
            i_pct = find(dt2<=3);% all bins with an ISI < 3ms
            
            pct = 0;%calculate pct isi
            if ~isempty(n_isi)
                [pct] = sum(n_isi(i_pct))/sum(n_isi);
                pct = round(pct*1000)/1000*100;
            end;
            
            %% check if cluster fulfills criteria
            if (pct < 3) && (mV > 40)
                
                %% make the PSTH
                n = zeros(length(dt),length(ts_trl));
                for jt = 1:length(ts_trl)
                    
                    n(:,jt) = hist(ts_trl{jt},dt);% count the number of spikes per time bin
                end;
                n(n>1)=1;
                
                
                %%                
                if ~isempty(n) && (median(sum(n,1))~=0)% if the PSTH is not empty
                    
                    %% save the channel info
                    ct = ct+1;a
                    CSC(ct,1) = spike_dat.label;
                    CSC(ct,2) = {gt};
                    
                    %perform interpolation of LFP in time domain to
                    %eliminate spike shape
                    int_step = [2 8];
                    for jt = 1:size(dum.trial,1)
                        lfp = squeeze( dum.trial(jt,:,:) );
                        spk = n(:,jt);
                        [lfp_interp] = interpLFP(lfp,spk,int_step);
                        
                        dum.trial(jt,:,:) = lfp_interp;
                    end;
                                        
                    %% SFC
                    i_base = find( dt >= -2000 & dt <0);
                    i_cue = find( dt >= 0 & dt <2000);
                    i_enc = find( dt >= 2000 & dt <4000);
                    
                    %extract LFP for different task intervals
                    lfp = squeeze(dum.trial)';
                    lfp_base = squeeze(dum.trial(:,1,i_base))';
                    lfp_cue = squeeze(dum.trial(:,1,i_cue))';
                    lfp_enc = squeeze(dum.trial(:,1,i_enc))';
                    
                    TW = 12;
                    k = 2*TW-1;
                    
                    params                  = [];
                    params.Fs               = 1e3;
                    params.pad              = -1;
                    params.fpass            = [1 100];
                    params.tapers           = [TW k];
                    params.err              = [2 0.05];
                    params.trialave         = 1;
                    
                    [~,~,~,~,~,f2,~,~,~,dum1]=coherencycpb(lfp_base,n(i_base,:),params,0);
                    [~,~,~,~,~,f2,~,~,~,dum2]=coherencycpb(lfp_cue,n(i_cue,:),params,0);
                    [~,~,~,~,~,f2,~,~,~,dum3]=coherencycpb(lfp_enc,n(i_enc,:),params,0);
                    
                    chck = [];
                    chck(1) =  max(max(dum1)) > .2;
                    chck(2) =  max(max(dum2)) > .2;                    
                    chck(3) =  max(max(dum3)) > .2;
                    
                    if any(chck ==1)
                        %% specgram of the SFC
                        movingwin = [.5 .025];
                        
                        TW = 2;
                        k = 2*TW-1;
                        
                        params                  = [];
                        params.Fs               = 1e3;
                        params.pad              = -1;
                        params.fpass            = [0 100];
                        params.tapers           = [TW k];
                        params.err              = [2 0.05];
                        params.trialave         = 1;
                        
                        y = squeeze(dum.trial)';
                        [Cgrm{ct},~,~,~,~,tAx,fAx]=cohgramcpb(y,n,movingwin,params,0);
                    end;
                end;
            end;
        end;
    end;
    
    fprintf('\n');
end;
