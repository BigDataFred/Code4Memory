
%% set Matlab path
restoredefaultpath;
addpath(genpath('/media/rouxf/rds-share/Fred/code/mcode/custom/'));
addpath(genpath('/media/rouxf/rds-share/Fred/code/mcode/toolboxes/chronux_2_11/'));
addpath(genpath('/media/rouxf/rds-share/Fred/code/mcode/toolboxes/osort-v3-rel/'));
addpath('/media/rouxf/rds-share/Fred/code/mcode/toolboxes/fieldtrip-20161009/');
ft_defaults;

%% reset workspace
clear all; close all; clc;

pID = 'P23AMS';% patient ID must be set manually
fprintf(['extracting cluster data for ',pID,'\n']);
fprintf([date,'\n']);

%% extract the source-path information and directory names
rpath = ['',filesep,'media',filesep,'rouxf',filesep,'rds-share',filesep,'iEEG_DATA',filesep,'MICRO',filesep,pID,filesep];

[exp,sesh] = get_EM_rec_info(rpath);

[rpath] = [rpath,exp,filesep];
[savepath] = rpath;

exp
sesh

%%
for gt = 1:length(sesh)%loop over the seshions
    
    %% load the Log data
    path_2_log = [rpath,sesh{gt},filesep,'log_dat',filesep];
    fn = dir([path_2_log,'*_EMtask_LogDat.mat']);
    
    [LogDat] = load([path_2_log,fn.name]);%load the data
    
    %% extract the trial information for the encoding phase
    [trl_all] = extract_Enc_trl(LogDat.LogDat1);% all, remembered (r), not-remembered (nr)
    
    %% initialize
    ct = 0;%counter
    clusterInfo =[];%save labels of significant channels

    %% spk-data source path and data files
    path_2_spk = [rpath,sesh{gt},filesep,'spike_dat',filesep];
    spk_files = dir([path_2_spk,'*.mat']);
    
    %%   load the spike data for the selected seshion
    for zt = 1:length(spk_files)
        
        fprintf([num2str(zt),filesep,num2str(length(spk_files))]);
        
        % load the spike data
        [dat] = load([path_2_spk,spk_files(zt).name]);
        
        [spike_dat] = dat.save_data{1}{1}{1};% structure with spike data
        fprintf('\n');
        
        clear dat;
        
        %%
        for it = 1:length(spike_dat.label) % loop over all clusters (ignore zero cluster)
            
            %% check if there are any spike times at all, if not ignore
            if ~isempty(spike_dat.timestamp{it}) && ( length(spike_dat.timestamp{it}) > 1 )
                
                %% increment counter
                ct = ct+1;
                
                %% extract spike times and compute pct of ISIs below 3ms               
                timestamps = spike_dat.timestamp{it};
                
                isi = diff(timestamps);
                isi = isi./1e3;% convert from micro-s to ms
                pctBelow = length(find(isi<=3.0))*100/length(isi);%pct below 3 ms
                
                %% compute maximal Voltage and SNR
                snr = mean(sqrt(median(spike_dat.waveform{it},3).^2))./mean(spike_dat.std(:,1));                
                
                %% compute power spectrum of spike times and spike-time autocorrelation
                spike_stats = getStatsForCluster([],timestamps);
                
                %% check power spectrum and auto-correlation
                [isOkPxxn] = checkPowerspectrum(spike_stats.Pxxn,spike_stats.f, 20, 100);%check for line noise
                
                sel_ix =(length(spike_stats.tvect)-1)/2+2:length(spike_stats.tvect);
                nnz =length(find(diff(round(spike_stats.Cxx(sel_ix).*1e3)/1e3)==0))/length(sel_ix);
                [isOkCxx] = nnz < 0.95;
                
                %% save the cluster stats and channel info                
                clusterInfo(ct,1:4) = [gt zt it-1 ct];% sesh microwire cluster counter
                
                %% check for single unit characteristics
                if (pctBelow < 3) && (snr>1) && (isOkPxxn ) % && isOkCxx                    
                    % save the single Unit info
                    clusterInfo(ct,5) = 1;% probably a SU
                else
                    clusterInfo(ct,5) = 0;% probably not a SU
                end;                                
                
                %% organize spike times into trials
                [ts_trl]  = cell(1,length(trl_all));
                
                c = 0;
                for jt = 1:length(trl_all)%loop over the trials
                    i_trl = find(spike_dat.trial{it} == trl_all(jt));% find all the timestamps for selected trial
                    ts_trl{jt} = spike_dat.time{it}( i_trl ).*1e3;% save the timestamps in ms
                end;
                
                %% estimate firing rates
                dt = -1000:100:4000;
                psth = [];
                for ft = 1:length( ts_trl )
                    psth(ft,:) = histc(ts_trl{ft},dt);
                end;
                %psth = psth./(100/1000);
                
                X(:,1)  = sum(psth(:,find(dt <=0)),2);
                X(:,2)  = sum(psth(:,find(dt >=100)),2);
                
                [~,~,~,stat] = ttest(X(:,2),X(:,1));
                empT = stat.tstat;
                
                randT = zeros(1,2e3);
                for jt = 1:2e3
                    randX = X(randperm(length(X(:))));%random partition
                    randX = reshape(randX,[size(X,1) 2]);% organize random partition into 2 sets
                    [~,~,~,stat] = ttest(randX(:,2),randX(:,1));
                    randT(jt) = stat.tstat;
                end;
                pval = length(find(abs(randT) >= abs(empT)))/length(randT);
                
                if pval < 0.05
                    if sign(empT) ==1
                        clusterInfo(ct,6:7) = [1 pval];
                    else
                        clusterInfo(ct,6:7) = [-1 pval];
                    end;
                else
                    clusterInfo(ct,6:7) = [0 pval];
                end;
                
%                 %% do spike count for different time windows
%                 tint = [-1e3 0; 0 1e3; 1e3 2e3; 2e3 3e3; 3e3 4e3];%these are the time bins
%                 [sc] = zeros(length(ts_trl),size(tint,1));
%                 for ft = 1:length( ts_trl )
%                     
%                     x = ts_trl{ft};% spike times in ms
%                     
%                     for kt = 1:size(tint,1)
%                         ix = find(x >= tint(kt,1) & x < tint(kt,2));% ts occurring in time period
%                         sc(ft,kt) = length(ix);% spike-count for time window
%                     end;
%                 end;
                
%                 %% do significance checks
%                 % calculate a threshold based on mean and std of spike count across trials during baseline
%                 thr = mean(sc(:,1))+2*std(sc(:,1),0,1);
%                 
%                 [chck] = zeros(size(sc,2)-1,2);
%                 c=0;
%                 for nt = 2:size(sc,2)
%                     c = c+1;
%                     [chck(c,1),~] = ttest(sc(:,nt),sc(:,1),'alpha',1e-3,'tail','right');
%                     [chck(c,2),~] = ttest(sc(:,nt),sc(:,1),'alpha',1e-3,'tail','left');
%                     %chck(c,3) = median(sc(:,nt)) >= 2;
%                     %chck(c,4) = median(sc(:,nt)) > thr;
%                 end;
                
%                 %% check if cluster fulfills criteria
%                 if any(sum(chck(:,1:2),2)>=1)
%                     if any(chck(:,1)) && ~any(chck(:,2))
%                         clusterInfo(ct,6) = 1;% significant activation
%                     elseif ~any(chck(:,1)) && any(chck(:,2));
%                         clusterInfo(ct,6) = -1;% significant deactivation
%                     else
%                         clusterInfo(ct,6) =2;% other
%                     end;
%                 else
%                     clusterInfo(ct,6) = 0;% no significant activation/deactivation
%                 end;
                
                
            end;
        end;
        fprintf('\n');
    end;
    
    %%
    out_name = [savepath,'Cluster_info_EM_',pID,'_',sesh{gt},'.mat'];
    readme = {'sessionID','microwireID','clusterID','unitID','spikeCount'};
    save(out_name,'clusterInfo','readme');
    
end;




