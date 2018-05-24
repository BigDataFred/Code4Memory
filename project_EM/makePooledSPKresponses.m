%% include in path
addpath(genpath('/media/rouxf/rds-share/Fred/code/mcode/custom/'));
addpath(genpath('/media/rouxf/rds-share/Fred/code/mcode/toolboxes/wave_clus-testing/'));

%% IDs of patients to include
pID = {'P02' 'P04' 'P05' 'P22AMS' 'P23AMS'};%patient ID

%% loop over patients
for pt = 1:length(pID)
    
    fprintf([num2str(pt),'/',num2str( length(pID) )],'\n');
    
    %% set the path for data access
    rpath = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',pID{pt},'/'];%rooth path
    
    [exp,sesh] = get_EM_rec_info(rpath);% extract session labels
    
    [rpath] = [rpath,exp,filesep]; % path to read data
    [savepath] = rpath;% path to save data
    
    %% initialize some variables for later use
    cluRaster = {};% spike raster
    FR = {};% firing rate
    spikes  = {};% waveforms
    classes = {};% cluster class
    clusterInfo = {};% cluster info
    sel_a = {};% selection index for behnke fried microwire
    sel_b = {};% selection index for cluster
    selIdxU = {};% selection index for clusters
    LogDat = {};
    trl_all = {};
    FRce = {};
    FRc = {};
    FRe = {};
    selIdxCE = {};
    selIdxC = {};
    selIdxE = {};
    nItemsCE = [];
    nItemsC = [];
    nItemsE = [];
    c1 = 0;
    c2 = 0;
    c3 = 0;
    
    %% extract responsive clusters and pool their firing rates
    for xt = 1:length( sesh )% loop across sessions
        
        %% load the spike-data
        path_2_spk = [rpath,sesh{xt},'/spike_dat/'];% path 2 spike data
        
        spk_files = dir([path_2_spk , '*spike_data*_stimlocked.mat' ]);% find files with spike data
        
        spike_dat = cell(1, length( spk_files ));
        for zt = 1:length( spk_files )%loop over the files with spike data
            
            [dat] = load([path_2_spk,spk_files(zt).name]);% load spike data
            [spike_dat{zt}] = dat.save_data{1}{1}{1};% extract fieldtrip-like structure with spike data
            
        end;
        clear dat;
        
        %% load the log-file data
        path_2_log = [rpath,sesh{xt},'/log_dat/'];% path to log-file data
        
        log_file = dir([path_2_log ,'*_Log*ile_EMtask_LogDat.mat']);
        
        LogDat{xt} = load([ path_2_log , log_file.name ]);
        
        [trl_all{xt}] = extract_Enc_trl(LogDat{xt}.LogDat1);
        nLog = length( trl_all{xt} );
        
        %% load the cluster-info
        p2d = rpath;
        dat = load([p2d,'Cluster_info_EM_',pID{pt},'_',sesh{xt},'.mat']);
        clusterInfo{xt} = dat.clusterInfo;
        CLUlabel{xt} = dat.CLUlabel;
        CSClabel{xt} = dat.CSClabel;
        
        %% check for changes in cluster classification
        p2cI = rpath;
        
        for it = 1:length( CLUlabel{xt} )
            
            fn = [CSClabel{xt}{it},'_',CLUlabel{xt}{it},'_fig1.tif'];
            
            [isSU] = dir([p2cI,sesh{xt},filesep,'cluster_figs/SU/',fn]);
            [isSUMU] = dir([p2cI,sesh{xt},filesep,'cluster_figs/SUMU/',fn]);
            [isMU] = dir([p2cI,sesh{xt},filesep,'cluster_figs/MU/',fn]);
            [isNoise] = dir([p2cI,sesh{xt},filesep,'cluster_figs/noise/',fn]);
            
            if ~isempty(isSU)
                clusterInfo{xt}(it,6) = 3;
            elseif ~isempty(isSUMU)
                clusterInfo{xt}(it,6) = 2;
            elseif ~isempty(isMU)
                clusterInfo{xt}(it,6) = 1;
            elseif  ~isempty(isNoise)
                clusterInfo{xt}(it,6) = 0;
            end;
            
        end;
        
        %% selection indexes for microwires and clusters
        selIdx1 = [find(clusterInfo{xt}(:,5)==1);find(clusterInfo{xt}(:,5)==2);find(clusterInfo{xt}(:,5)==3)];% get all clusters with SU, MU and SUMU characteristics
        selIdx2 = [find(clusterInfo{xt}(:,12)==1);find(clusterInfo{xt}(:,14)==1)];% get all clusters that show significant change in firing rate
        
        selIdxU{xt} = intersect( selIdx1 , selIdx2 );% keep intersection of both selection indexes
        
        %% make a spike-raster for each selected cluster
        sel_a{xt} = clusterInfo{xt}(selIdxU{xt},2);% microwire index
        sel_b{xt} = clusterInfo{xt}(selIdxU{xt},3);% cluster index
        
        %%
        if ~isempty( sel_a{xt} )
            % if ~isempty( sel_a{xt} )
            %     chck = [sel_a{xt} sel_b{xt}];
            %     idx = [];
            %     for it = 1:size(chck,1)
            %         idx{it} = find(ismember(chck,chck(it,:),'rows'));
            %         idx{it} = idx{it}(1);
            %     end;
            %
            %     idx = unique([idx{:}]);
            %
            %     chck = chck(idx,:);
            %
            %     sel_a{xt} = chck(:,1);
            %     sel_b{xt} = chck(:,2);
            % end;
            
            dt = -2000:1:5000;% time bins of spike raster ( should have 1ms width )
            Fs = 1/((dt(2)-dt(1))/1000);% sampling frequency
            
            for it = 1:length( sel_a{xt} )
                
                ts = spike_dat{sel_a{xt}(it)}.time{sel_b{xt}(it)+1};%extract timestamps
                ts = ts.*1e3;% convert spike-times from s into ms
                trl = spike_dat{sel_a{xt}(it)}.trial{sel_b{xt}(it)+1};% these are the trial labels for each spike time
                
                if ~isempty(ts)
                    x=[];
                    for jt = 1:length(trl_all{xt})%loop across trials
                        
                        ix = find(trl == trl_all{xt}(jt));
                        x = ts(ix);% readout the relevant spike-times
                        
                        x(x<=dt(1)) = [];
                        x(x>dt(end)) = [];
                        
                        cluRaster{xt}{it}(jt,:) = hist(x,dt);% 0= no spike, 1 = spike
                        
                    end;
                else
                    cluRaster{xt}{it} = [];
                end;
                
            end;
            
            %% compute firing rate
            T = 0.5; % this is the width over which firing rate should be integrated (in s)
            step = T*Fs;% convert T into samples/bins
            
            for ct = 1:length( cluRaster{xt} )
                
                tInt = [-1000:step:0];% baseline time window (work in ms)
                x = [];
                nSpk = 0;
                for ft = 1:length(tInt)-1
                    nSpk = sum(cluRaster{xt}{ct}(:,find(dt >=tInt(ft) & dt < tInt(ft+1))),2);
                    x(:,ft) = nSpk./T;% firing rate
                end;
                FR{xt}{ct}(:,1) = mean(x,2);
                clear x;
                
                tInt = [0:step:2000];% cue time windoe
                x = [];
                nSpk = 0;
                for ft = 1:length(tInt)-1
                    nSpk = sum(cluRaster{xt}{ct}(:,find(dt >=tInt(ft) & dt < tInt(ft+1))),2);
                    x(:,ft) = nSpk./T;% firing rate
                end;
                FR{xt}{ct}(:,2) = mean(x,2);
                clear x;
                
                tInt = [2000:step:5000];% encoding time window
                x = [];
                nSpk = 0;
                for ft = 1:length(tInt)-1
                    nSpk = sum(cluRaster{xt}{ct}(:,find(dt >=tInt(ft) & dt < tInt(ft+1))),2);
                    x(:,ft) = nSpk./T;%firing rate
                end;
                FR{xt}{ct}(:,3) = mean(x,2);
                clear x;
                
            end;
            
            %% save the waveform and class-labels of each cluster
            for  jt = 1:length( sel_a{xt} )
                
                spikes{xt}{jt} = [];
                classes{xt}{jt} = [];
                
                for it = 1:length(spike_dat{sel_a{xt}(jt)}.waveform)
                    wvf = squeeze(spike_dat{sel_a{xt}(jt)}.waveform{it})';
                    if size(wvf,2) ==1
                        wvf = wvf';
                    end;
                    
                    spikes{xt}{jt} = [spikes{xt}{jt};wvf];
                    classes{xt}{jt} = [classes{xt}{jt};it*ones(size(wvf,1),1)];
                    
                end;
            end;
            
            %% SME          
            selIdxCE{xt} = intersect(find(clusterInfo{xt}(selIdxU{xt},12)),find(clusterInfo{xt}(selIdxU{xt},14)));
            selIdxC{xt} = setdiff(find(clusterInfo{xt}(selIdxU{xt},12)), selIdxCE{xt});
            selIdxE{xt} = setdiff(find(clusterInfo{xt}(selIdxU{xt},14)), selIdxCE{xt});
            
            n = [];
            n(1) = length(LogDat{xt}.ix{4});
            n(2) = length(LogDat{xt}.ix{5});
            n(3) = length(LogDat{xt}.ix{6});
            
            if sum(n) ~= length( trl_all{xt} )
                error('trial count must match');
            end;
            
            cond = zeros(1,length( trl_all{xt}) );
            cond( LogDat{xt}.ix{4} ) = 2;
            cond( LogDat{xt}.ix{5} ) = 1;
            cond( LogDat{xt}.ix{6} ) = 0;
            
            for it = 1:length(selIdxCE{xt})
                c1 = c1+1;
                
                FRce{c1} = mean(FR{xt}{selIdxCE{xt}(it)}(:,[1 2 3]),1);% pool the firing rates across all sessions and clusters
                
                ix1 = find(mean(FR{xt}{selIdxCE{xt}(it)}(:,2:3),2) > median(mean(FR{xt}{selIdxCE{xt}(it)}(:,2:3)),2));
                ix2 = find(mean(FR{xt}{selIdxCE{xt}(it)}(:,2:3),2) <= median(mean(FR{xt}{selIdxCE{xt}(it)}(:,2:3)),2));
                
                nItemsCE(c1,1) = mean(cond(ix2));
                nItemsCE(c1,2) = mean(cond(ix1));
            end;
            
            for it = 1:length(selIdxC{xt})
                c2 = c2+1;
                
                FRc{c2} = mean(FR{xt}{selIdxC{xt}(it)}(:,[1 2]),1);% pool the firing rates across all sessions and clusters
                
                ix1 = find(FR{xt}{selIdxC{xt}(it)}(:,2) > median(FR{xt}{selIdxC{xt}(it)}(:,2)));
                ix2 = find(FR{xt}{selIdxC{xt}(it)}(:,2) <= median(FR{xt}{selIdxC{xt}(it)}(:,2)));
                
                nItemsC(c2,1) = mean(cond(ix2));
                nItemsC(c2,2) = mean(cond(ix1));
            end;
            
            for it = 1:length(selIdxE{xt})
                c3 = c3+1;
                
                FRe{c3} = mean(FR{xt}{selIdxE{xt}(it)}(:,[1 3]),1);% pool the firing rates across all sessions and clusters
                
                ix1 = find(FR{xt}{selIdxE{xt}(it)}(:,3) > median(FR{xt}{selIdxE{xt}(it)}(:,3)));
                ix2 = find(FR{xt}{selIdxE{xt}(it)}(:,3) <= median(FR{xt}{selIdxE{xt}(it)}(:,3)));
                
                nItemsE(c3,1) = mean(cond(ix2));
                nItemsE(c3,2) = mean(cond(ix1));
            end;
            
        end;
    end;
    
    %% save the pooled data
    dat = [];
    dat.cluRaster = cluRaster;% spike raster
    dat.FR = FR;% firing rate
    dat.spikes  = spikes;% waveforms
    dat.classes = classes;% cluster class
    dat.clusterInfo = clusterInfo;% cluster info
    dat.sel_a = sel_a;% selection index for behnke fried microwire
    dat.sel_b = sel_b;% selection index for cluster
    dat.selIdxU = selIdxU;% selection index for clusters
    dat.LogDat = LogDat;
    dat.trl_all = trl_all;
    dat.dt = dt;
    dat.Fs = Fs;
    dat.FRce = FRce;
    dat.FRc = FRc;
    dat.FRe = FRe;
    dat.nItemsCE = nItemsCE;
    dat.nItemsC = nItemsC;
    dat.nItemsE = nItemsE;
    dat.selIdxCE = selIdxCE;
    dat.selIdxC = selIdxC;
    dat.selIdxE = selIdxE;
    dat.CSC = CSClabel;
    
    save([savepath,pID{pt},'_pooledSPKdataEMtask.mat'],'dat');
    
end;

%%
clear all;
exit;
