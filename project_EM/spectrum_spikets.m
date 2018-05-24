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
pID = 'P23AMS';% patient ID must be set manually
rpath = ['',filesep,'media',filesep,'rouxf',filesep,'rds-share',filesep,'iEEG_DATA',filesep,'MICRO',filesep,pID,filesep];

[exp,sesh] = get_EM_rec_info( rpath );
exp
sesh

rpath = [rpath,exp,filesep]
%%
ct = 0;%counter
CSC ={};%save labels of significant channels

for gt = 1:length(sesh)%loop over the seshions
    
    %% spk-data source path and data files
    path_2_spk = [rpath,sesh{gt},filesep,'spike_dat',filesep];
    spk_files = dir([path_2_spk,'*.mat']);
    
    %% load the Log data
    path_2_log = [rpath,sesh{gt},filesep,'log_dat',filesep];
    fn = dir([path_2_log,'*_LogFile_EMtask_LogDat.mat']);
    
    LogDat = load([path_2_log,fn.name]);%load the data
        
    %% extract the trial information for the encoding phase
    [trl_all,trl_r,trl_nr] = extract_Enc_trl(LogDat.LogDat1,LogDat.ix);% all, remembered (r), not-remembered (nr)
    
    %% load the LFP data for the selected seshion
    path_2_lfps = [rpath,sesh{gt},filesep,'lfp_dat',filesep];
    lfp_files = dir([path_2_lfps,'*.mat']);
    
    [lfp_dat] = load_lfp_data(lfp_files,path_2_lfps,trl_all);
    
    %%    loop over all CSC channels
    for zt = 1:length(spk_files)
        
        fprintf([num2str(zt),filesep,num2str(length(spk_files))]);
        
        % load the spike data
        [dat] = load([path_2_spk,spk_files(zt).name]);
        
        [spike_dat] = dat.save_data{1}{1}{1};% structure with spike data
        clear dat;
        
        %%
        dt = -2000:1:4000;% used for the binning of spike times (PSTH)   
        
        i_base = find( dt >= -2e3 & dt < 0 );
        i_cue = find( dt >= 0 & dt < 2e3 );
        i_enc = find( dt >= 2e3 & dt < 4e3 );
        
        Fs = 1e3; % Sampling Rate of the binned point process
        
        cID = unique(spike_dat.unit{1});% IDs of individual clusters
        
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
            
            dt2 = 0:1.001e3;%bins in ms
            [n_isi] = hist(isi,dt2);
            n_isi(end) = [];
            dt2(end) = [];
            
            i_pct = find(dt2<=3);% all bins with an ISI < 3ms
            
            pct = 0;%calculate pct isi
            if ~isempty(n_isi)
                [pct] = sum(n_isi(i_pct))/sum(n_isi);
                pct = round(pct*1000)/1000*100;
            end;                                                      
            
            %% search for significant difference between the number of spikes
            sc1 = zeros(length(ts_trl),1);
            sc2 = zeros(length(ts_trl),1);
            sc3 = zeros(length(ts_trl),1);
            
            for ft = 1:length( ts_trl )
                   
                x = ts_trl{ft};
                
                n_base = find(x >= -1e3 & x < 0);% ts in baseline period
                n_enc1 = find(x >= 0 & x < 1e3); % ts in 1st half of encoding
                n_enc2 = find(x >= 2e3 & x < 3e3);% ts in 2nd half of encoding
                
                sc1(ft) = length(n_base);% spike-count baseline
                sc2(ft) = length(n_enc1);% spike-count encoding 1
                sc3(ft) = length(n_enc2);% spike-count encoding 2
                
            end;
            
            %% do significance checks
            % calculate a threshold based on mean and std of spike count across trials during baseline
            thr = mean(sc1)+2*std(sc1,0,1);
            
            chck1 = []; chck2 = [];
            [chck1(1),p(1)] = ttest(sc2,sc1,'alpha',1e-2,'tail','right');
            [chck2(1),p(2)] = ttest(sc3,sc1,'alpha',1e-2,'tail','right');
            
            % median number of spikes needs to be at least 2 or more
            chck1(2) = median(sc2) >= 2;
            chck2(2) = median(sc3) >= 2;
            
            %chck1(3) = median(sc2) > thr;
            %chck2(3) = median(sc3) > thr;
            
            %% check if cluster fulfills criteria
            if ((sum(chck1) ==2) ||  (sum(chck2) ==2)) && (pct < 3) && (mV > 20)
                
                %% make the PSTH
                n = zeros(length(dt),length(ts_trl));
                for jt = 1:length(ts_trl)                                       
                    
                    n(:,jt) = hist(ts_trl{jt},dt);% count the number of spikes per time bin
                end;
                n(n>1)=1;
                
                [thr2] = mean(sum(n(find(dt >=-1e3 & dt <0),:),1))+2*std(sum(n(find(dt >=-1e3 & dt <0),:),1),0,2);
                if floor(thr*10)/10 ~= floor(thr2*10)/10
                    error('threshold estimates must be equal');
                end;
                
                %% smooth firing rate for visualization
                nw = 50;
                w = gausswin(21);%gaussian kernel
                dt3 = -2e3:nw:4e3;% time bins with larger width
                n2 = zeros(length(dt3),length(ts_trl));
                for jt = 1:length(ts_trl)
                    dum = hist(ts_trl{jt},dt3);
                    dum(dum>1)=1;
                    sn = length(dum);
                    dum = conv([fliplr(dum) dum fliplr(dum)],w,'same')./sum(w);%smooth the spike count over time
                    n2(:,jt) = dum(sn+1:2*sn);
                end;
                
                %%
                if ~isempty(n)% if the PSTH is not empty
                    
                    %% save the channel info
                    ct = ct+1;
                    CSC(ct,1) = spike_dat.label;
                    CSC(ct,2) = {gt};         
                    
                    %%
                    isi = diff([ts_trl{:}]);
                    isi = isi(isi>0);
                    
                    dt2 = 0:1001;%bins in ms
                    [n_isi] = hist(isi,dt2);
                    n_isi(end) = [];
                    dt2(end) = [];
            
                    %% autocorrelation for each task interval                   
                    
                    xcn = 1e3;
                    xc1 = zeros(xcn*2+1,size(n,2));
                    xc2 = zeros(xcn*2+1,size(n,2));
                    xc3 = zeros(xcn*2+1,size(n,2));
                    for ft = 1:size(n,2)
                        x = n(i_base,ft);%baseline
                        [xc1(:,ft),lag] = xcorr(x-mean(x),xcn,'coeff');
                        x = n(i_cue,ft);%cue
                        [xc2(:,ft),lag] = xcorr(x-mean(x),xcn,'coeff');
                        x = n(i_enc,ft);%encoding
                        [xc3(:,ft),lag] = xcorr(x-mean(x),xcn,'coeff');
                    end;                                        
                    
                    %% LFP data
                    cfg                     = [];
                    cfg.channel             = spike_dat.label;
                    
                    [dum] = ft_selectdata( cfg , lfp_dat );%select the relevant CSC channel
                    
                    %perform interpolation of LFP in time domain to
                    %eliminate spike shape
                    int_step = [2 8];
                    for jt = 1:size(dum.trial,1)
                        lfp = squeeze( dum.trial(jt,:,:) );
                        spk = n(:,jt);
                        [lfp_interp] = interpLFP(lfp,spk,int_step);
                        
                         dum.trial(jt,:,:) = lfp_interp;
                    end;
                    
                    %extract LFP for different task intervals
                    lfp_base = squeeze(dum.trial(:,1,i_base))';
                    lfp_cue = squeeze(dum.trial(:,1,i_cue))';
                    lfp_enc = squeeze(dum.trial(:,1,i_enc))';
                    
                    %% STA
                    win = 0.25e3;
                    
                    STA = zeros(win*2+1,size(n,2));
                    c = 0;
                    for jt = 1:size(n,2)%loop over trials
                        ix = find(n(:,jt)~=0);% find spike time indexes
                        for kt = 1:length(ix)%loop over number of spikes
                            if ((ix(kt))-win >0) && ((ix(kt))+win <=size(n,1))
                                x = squeeze(dum.trial(jt,:,ix(kt)-win:ix(kt)+win));
                                STA(:,jt) = STA(:,jt) + x;
                            end;
                        end;
                    end;
                    STA = STA./length(ix);
                    
%                     STA1 = zeros(win*2+1,size(n,2));
%                     c = 0;
%                     for jt = 1:size(n,2)%loop over trials
%                         ix = find(n(i_base,jt)~=0);% find spike time indexes
%                         for kt = 1:length(ix)%loop over number of spikes
%                             if ((ix(kt))-win >0) && ((ix(kt))+win <=size(n,1))
%                                 x = squeeze(dum.trial(jt,:,ix(kt)-win:ix(kt)+win));
%                                 STA1(:,jt) = STA1(:,jt) + x;
%                             end;
%                         end;
%                     end;
%                     STA1 = STA1./length(ix);
%                                                             
%                     STA2 = zeros(win*2+1,size(n,2));
%                     c = 0;
%                     for jt = 1:size(n,2)%loop over trials
%                         ix = find(n(i_cue,jt)~=0);% find spike time indexes
%                         for kt = 1:length(ix)%loop over number of spikes
%                             if ((ix(kt))-win >0) && ((ix(kt))+win <=size(n,1))
%                                 x = squeeze(dum.trial(jt,:,ix(kt)-win:ix(kt)+win));
%                                 STA2(:,jt) = STA2(:,jt) + x;
%                             end;
%                         end;
%                     end;
%                     STA2 = STA2./length(ix);
%                                                             
%                     STA3 = zeros(win*2+1,size(n,2));
%                     c = 0;
%                     for jt = 1:size(n,2)%loop over trials
%                         ix = find(n(i_enc,jt)~=0);% find spike time indexes
%                         for kt = 1:length(ix)%loop over number of spikes
%                             if ((ix(kt))-win >0) && ((ix(kt))+win <=size(n,1))
%                                 x = squeeze(dum.trial(jt,:,ix(kt)-win:ix(kt)+win));
%                                 STA3(:,jt) = STA3(:,jt) + x;
%                             end;
%                         end;
%                     end;
%                     STA3 = STA3./length(ix);
                    
                    %% spectrum and SFC
                    TW = 6;
                    k = 2*TW-1;
                    
                    params                  = [];
                    params.Fs               = 1e3;
                    params.pad              = -1;
                    params.fpass            = [1 100];
                    params.tapers           = [TW k];
                    params.err              = [2 0.05];
                    params.trialave         = 1;
                    
                    [Snn1,f1,R1,Serr1]=mtspectrumpb(n(i_base,:),params,0);
                    [Snn2,f1,R2,Serr2]=mtspectrumpb(n(i_cue,:),params,0);
                    [Snn3,f1,R3,Serr3]=mtspectrumpb(n(i_enc,:),params,0);
                    
                    %[~,dy] = gradient(lfp_base);
                    [Syy1,f1,Serr4]=mtspectrumc(lfp_base,params);
                    %[~,dy] = gradient(lfp_cue);
                    [Syy2,f1,Serr5]=mtspectrumc(lfp_cue,params);
                    %[~,dy] = gradient(lfp_enc);
                    [Syy3,f1,Serr6]=mtspectrumc(lfp_enc,params);
                    
                    [Syy1] = correct_one_over_f(Syy1, f1);
                    [Syy2] = correct_one_over_f(Syy2, f1);
                    [Syy3] = correct_one_over_f(Syy3, f1);
                    
                    TW = 12;
                    k = 2*TW-1;
                    
                    params                  = [];
                    params.Fs               = 1e3;
                    params.pad              = -1;
                    params.fpass            = [1 100];
                    params.tapers           = [TW k];
                    params.err              = [2 0.05];
                    params.trialave         = 1;
                    
                    [C1,~,~,~,~,f2]=coherencycpb(lfp_base,n(i_base,:),params,0);
                    [C2,~,~,~,~,f2]=coherencycpb(lfp_cue,n(i_cue,:),params,0);
                    [C3,~,~,~,~,f2]=coherencycpb(lfp_enc,n(i_enc,:),params,0);
                    
                    %% specgram of the STA
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
                    [Cgrm,~,~,~,~,tAx,fAx]=cohgramcpb(y,n,movingwin,params,0);
                    
                    
                    
                end;
            end;
        end;
        
        fprintf('\n');
    end;
end;

