%% reset workspace
clear all; close all; clc;

pID = 'P23AMS';% patient ID must be set manually
fprintf(['running spectrum_spikets for ',pID,'\n']);
fprintf([date,'\n']);

%% set Matlab path
restoredefaultpath;
addpath(genpath('/media/rouxf/rds-share/Fred/code/mcode/custom/'));
addpath(genpath('/media/rouxf/rds-share/Fred/code/mcode/toolboxes/chronux_2_11/'));
addpath(genpath('/media/rouxf/rds-share/Fred/code/mcode/toolboxes/osort-v3-rel/'));
addpath('/media/rouxf/rds-share/Fred/code/mcode/toolboxes/fieldtrip-20161009/');
ft_defaults;

%% extract the source-path information and directory names
rpath = ['',filesep,'media',filesep,'rouxf',filesep,'rds-share',filesep,'iEEG_DATA',filesep,'MICRO',filesep,pID,filesep];
savepath = rpath;

[exp,sesh] = get_EM_rec_info(rpath);
exp
sesh
rpath = [rpath,exp,filesep]

%% initialize
ct = 0;%counter
CSC =[];%save labels of significant channels
%%
for gt = 1:length(sesh)%loop over the seshions
    
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
        for it = 2:length(spike_dat.label) % loop over all clusters (ignore zero cluster)
            
            if ~isempty(spike_dat.timestamp{it}) && ( length(spike_dat.timestamp{it}) > 1 )
                %% compute maximal Voltage and SNR
                
                timestamps = spike_dat.timestamp{it};
                
                [spike_stats] = getStatsForCluster([], timestamps);
                
                mV = max(abs(mean(spike_dat.waveform{it},3)));
                snr = mean(sqrt(median(spike_dat.waveform{it},3).^2))./mean(spike_dat.std(:,1));
                nevents = size(squeeze(spike_dat.waveform{it}),2);
                
                %% check power spectrum and auto-correlation
                [isOkPxxn] = checkPowerspectrum(spike_stats.Pxxn,spike_stats.f, 20, 500);%check for line noise
                [isOkCxx] = length(find(round(diff(spike_stats.Cxx)*1e6)/1e6==0)) ==0;
                
                %% check for single unit characteristics
                if (spike_stats.percentageBelow(3) < 3) && (snr>1) && (isOkPxxn && isOkCxx)
                    
                    %% save the cluster stats and channel info
                    ct = ct+1;
                    CSC(ct,1:4) = [gt zt it ct];% sesh microwire cluster counter
                    
                    %% load the Log data
                    path_2_log = [rpath,sesh{gt},filesep,'log_dat',filesep];
                    fn = dir([path_2_log,'*_EMtask_LogDat.mat']);
                    
                    LogDat = load([path_2_log,fn.name]);%load the data
                    
                    %% extract the trial information for the encoding phase
                    [trl_all] = extract_Enc_trl(LogDat.LogDat1,LogDat.ix);% all, remembered (r), not-remembered (nr)
                    
                    %% organize spike times into trials
                    ts_trl  = cell(1,length(trl_all));
                    
                    c = 0;
                    for jt = 1:length(trl_all)%loop over the trials
                        i_trl = find(spike_dat.trial{it} == trl_all(jt));% find all the timestamps for selected trial
                        ts_trl{jt} = spike_dat.time{it}( i_trl ).*1e3;% save the timestamps in ms
                    end;
                    
                    %% search for significant difference between the number of spikes
                    sc = zeros(length(ts_trl),5);
                    tint = [-1e3 0; 0 1e3; 1e3 2e3; 2e3 3e3; 3e3 4e3];
                    for ft = 1:length( ts_trl )
                        
                        x = ts_trl{ft};% spike times in ms
                        
                        for kt = 1:size(tint,1)
                            ix = find(x >= tint(kt,1) & x < tint(kt,2));% ts in baseline period
                            sc(ft,kt) = length(ix);% spike-count baseline
                        end;
                    end;
                    
                    %% do significance checks
                    % calculate a threshold based on mean and std of spike count across trials during baseline
                    thr = mean(sc(:,1))+2*std(sc(:,1),0,1);
                    
                    chck = zeros(size(sc,2)-1,3);
                    c=0;
                    for nt = 2:size(sc,2)
                        c = c+1;
                        [chck(c,1),~] = ttest(sc(:,nt),sc(:,1),'alpha',1e-2,'tail','right');
                        chck(c,2) = median(sc(:,nt)) >= 2;
                        chck(c,3) = median(sc(:,nt)) > thr;
                    end;
                    
                    %% check if cluster fulfills criteria
                    if any(sum(chck,2)>=2)
                        CSC(ct,5) = 1;% significant activation
                    else
                        CSC(ct,5) = 0;% no activation
                    end;
                    
                end;
            end;
        end;
        fprintf('\n');
    end;
end;

%%
out_name = ['Cluster_info_',pID,',mat'];
save([savepath,out_name],'CSC');


