%% set Matlab path
restoredefaultpath;
addpath(genpath('/home/rouxf/prj/Bham/code/mcode/'));
addpath(genpath('/home/rouxf/tbx/chronux_2_11/'));
addpath(genpath('/home/rouxf/tbx/osort-v3-rel/'));

%%
if isempty(gcp('NoCreate'))
    parpool(36,'SpmdEnabled',false);
end;

%% reset workspace
bw = 0.5;%set bin width for firing rate
nperms = 2e3;% permutation for stats

pID = {'P22AMS','P23AMS','P02','P04','P05','P07','P08'};% patient ID must be set manually

%%
for st = 6:length(pID)
    
    fprintf(['extracting cluster data for ',pID{st},'\n']);
    fprintf([date,'\n']);
    
    %% extract the source-path information and directory names
    rpath = ['',filesep,'media',filesep,'rouxf',filesep,'rds-share',filesep,'iEEG_DATA',filesep,'MICRO',filesep,pID{st},filesep];
    
    [exp,sesh] = get_EM_rec_info(rpath);% get the experiment foldername and foldername of each session
    
    %%
    for et = 1:length(exp)
        
        [rpath2] = [rpath,exp{et},filesep];
        [savepath] = rpath2;
        
        exp{et}
        sesh{et}
        
        %%
        for gt = 1:length(sesh{et})%loop over the seshions
            
            [sesh2] = sesh{et}{gt};
            
            chck = dir([savepath,sesh2 ,filesep,'cluster_figs',filesep]);
            if isempty(chck)
                mkdir([savepath,sesh2 ,filesep,'cluster_figs',filesep]);
                mkdir([savepath,sesh2 ,filesep,'cluster_figs',filesep,'SU',filesep]);
                mkdir([savepath,sesh2 ,filesep,'cluster_figs',filesep,'SUMU',filesep]);
                mkdir([savepath,sesh2 ,filesep,'cluster_figs',filesep,'MU',filesep]);
                mkdir([savepath,sesh2 ,filesep,'cluster_figs',filesep,'noise',filesep]);
            else
                rmdir([savepath,sesh2 ,filesep,'cluster_figs',filesep],'s');
                mkdir([savepath,sesh2 ,filesep,'cluster_figs',filesep]);
                mkdir([savepath,sesh2 ,filesep,'cluster_figs',filesep,'SU',filesep]);
                mkdir([savepath,sesh2 ,filesep,'cluster_figs',filesep,'SUMU',filesep]);
                mkdir([savepath,sesh2 ,filesep,'cluster_figs',filesep,'MU',filesep]);
                mkdir([savepath,sesh2 ,filesep,'cluster_figs',filesep,'noise',filesep]);
            end;
            
            %% load the Log data
            try
                path_2_log = [savepath,sesh2 ,filesep,'log_dat',filesep];
                fn = dir([path_2_log,'*_EMtask_LogDat.mat']);
                [LogDat] = load([path_2_log,fn.name]);%load the data
                
                %% do some processing on the Log data
                for it = 1:length(LogDat.ix)
                    if size(LogDat.ix{it},2) > size(LogDat.ix{it},1)
                        LogDat.ix{it} = LogDat.ix{it}';
                    end;
                end;
                
                trlIx1 = [LogDat.ix{1};LogDat.ix{2};LogDat.ix{3}];% f-f, p-p, f-p
                trlIx2 = [LogDat.ix{4};LogDat.ix{5};LogDat.ix{6}];% hit, 1/2 , miss
                
                %count trials per condition
                nTrl = zeros(1,length( LogDat.ix ));
                for it = 1:length(LogDat.ix)
                    nTrl(it) = length(LogDat.ix{it});
                end;
                
                %% extract the trial information for the encoding phase
                [trl_all] = extract_Enc_trl(LogDat.LogDat1);% all, remembered (r), not-remembered (nr)
                
                %do some sanity checks on the trial-count
                if (sum(nTrl(1:3)) ~= length(trl_all)) || (sum(nTrl(4:6)) ~= length(trl_all))
                    error('number of trial events is out of range');
                end;
                
                if (length(trlIx1) ~= length(trl_all)) || (length(trlIx2) ~= length(trl_all))
                    error('number of trial events is out of range');
                end;
                
                %% initialize
                ct = 0;%counter incremented with each cluster
                clusterInfo =[];%save labels of significant channels
                FR = {};% firing rates
                empT1 = [];    empT2 = []; % stats
                empF1 = [];    empF2 = [];
                CSClabel = {};% cluster - mircowire lookup
                CLUlabel = {};% cluser name
                
                %% spk-data source path and data files
                path_2_spk = [savepath,sesh2 ,filesep,'spike_dat',filesep];
                spk_files = dir([path_2_spk,'*.mat']);
                
                %%   load the spike data for the selected session
                for zt = 1:length(spk_files)
                    
                    fprintf([num2str(zt),filesep,num2str(length(spk_files))]);
                    % load the spike data
                    load([path_2_spk,spk_files(zt).name]);
                    [spike_dat] = save_data{1}{1}{1};% structure with spike data
                    clear save_data;
                    fprintf('\n');
                    
                    %%
                    for it = 1:length(spike_dat.label) % loop over all clusters (including zero cluster)
                        
                        %% get label of micro-wire
                        rix1 = regexp(spk_files(zt).name,'spike_data_')+length('spike_data_');
                        rix2 = regexp(spk_files(zt).name,'\d{4}')-2;
                        chan = spk_files(zt).name( rix1:rix2 );
                        cluname = spike_dat.label{it};
                        
                        %% check if there are any spike times at all, if not ignore
                        if ~isempty(spike_dat.timestamp{it}) && ( length(spike_dat.timestamp{it}) > 1 )
                            
                            %% save the cluster stats and channel info
                            ct = ct+1;% increment counter
                            clusterInfo(ct,1:4) = [gt zt it-1 ct];% sesh microwire cluster counter
                            CSClabel = [CSClabel;chan];% microwire label
                            CLUlabel = [CLUlabel;cluname];
                            
                            %% extract spike times
                            timestamps = spike_dat.timestamp{it};
                            
                            %% spike train autocorrelation
                            x =timestamps./1e3; % convert micro-seconds to ms
                            dt = 0:1:max(x);
                            [n,~] = hist(x,dt);
                            
                            [Cxx,lag] = xcorr(n,1000,'coeff');
                            
                            % check the values of the spike-trrai autocorrelation
                            pctZero = length(find(diff(round(Cxx*1e4)/1e4)==0))/length(lag);
                            if  pctZero<0.95
                                isOkCxx = true;
                            else
                                isOkCxx = false;
                            end;
                            
                            %% spike train powerspectrum
                            [spike_stats] = getStatsForCluster(timestamps);
                            [isOkPxxn] = checkPowerspectrum(spike_stats.Pxxn,spike_stats.f, 20, 200);
                            
                            %% compute pct of ISIs below 3ms
                            isi = diff(timestamps);
                            isi = isi./1e3;% convert from micro-s to ms
                            [pctBelow] = length(find(isi<=3.0))*100/length(isi);%pct below 3 ms
                            [pctAbove] = length(find(isi>=1000.0))*100/length(isi);%pct below 3 ms
                            
                            %% extract the waveforms
                            wvf = squeeze(spike_dat.waveform{it});
                            [~,ix1]=max(abs(wvf),[],1);
                            
                            h = zeros(size(wvf'));
                            for jt = 1:size(wvf,2)
                                h(jt,:) = hilbert(wvf(:,jt))';
                                h(jt,:) = h(jt,:)./abs(h(jt,:));
                            end;
                            itc = abs(mean(h,1));
                            clear h;
                            [~,ix2] = max(itc);
                            
                            [wvfStab] = (floor(abs(ix2-mean(ix1))) == 0);
                            
                            %% measure polarity of waveshape at maximum
                            chck = [];
                            for kt = 1:length( ix1 )
                                chck(kt) = sign(wvf(ix1(kt),kt));
                            end;
                            nneg = length(find(chck==-1));
                            npos = length(find(chck==1));
                            if npos > nneg
                                [pctSign] = nneg/length(chck);
                            else
                                [pctSign] = npos/length(chck);
                            end;
                            
                            %% compute maximal SNR
                            [snr] = max(sqrt(median(wvf,1).^2))./mean(spike_dat.std(:,1));
                            
                            %% check for single unit characteristics
                            clucat = [];
                            if (pctBelow < 3 && pctAbove < 90) && (snr>1) && (isOkPxxn ) && (isOkCxx ) && (max(itc) >.9) && (pctSign <0.01) && (wvfStab)
                                % save the single Unit info
                                clusterInfo(ct,5) = 3;% probably a SU
                                plotClusterFeatures;
                                clucat = 'SU';
                            elseif (pctBelow < 10 && pctAbove < 90) && (snr>1) && (isOkPxxn ) && (isOkCxx ) && (max(itc) >.75) && (pctSign <0.1)
                                clusterInfo(ct,5) = 2;% probably between a SU and a MU
                                plotClusterFeatures;
                                clucat = 'SUMU';
                            elseif (pctAbove < 90) && (snr>1) && (isOkPxxn ) && (isOkCxx ) %&& (max(itc) >.3)
                                clusterInfo(ct,5) = 1;% probably a MU
                                plotClusterFeatures;
                                clucat = 'MU';
                            else
                                clusterInfo(ct,5) = 0;% probably noise
                                plotClusterFeatures;
                                clucat = 'noise';
                            end;
                            
                            nFig = get(gcf,'Number');
                            for ft = 1:nFig
                                print(ft,'-dtiff','-r100',[savepath,sesh2 ,filesep,'cluster_figs/',clucat,'/',spike_dat.hdr.label,'_',spike_dat.label{it},'_fig',num2str(ft),'.tif']);
                            end;
                            close all;
                            
                            % this is to correct the initial classification
                            clusterInfo(ct,6) = 0;
                            
                            %save the isi pct and snr
                            clusterInfo(ct,7) = pctBelow;
                            clusterInfo(ct,8) = snr;
                            clusterInfo(ct,9) = max(itc);
                            clusterInfo(ct,10) = pctSign;
                            
                            %% estimate spike count
                            dt = -2000:1:5000;% time bins in ms
                            Fs = 1/((dt(2)-dt(1))/1000);% set the sampling rate of the spike-raster matrix
                            
                            spkRaster = zeros(length(trl_all),length(dt));
                            for ft = 1:length(trl_all)%loop over trials
                                trlIx = find(spike_dat.trial{it} == trl_all(ft));%find indexes of trials
                                x = spike_dat.time{it}( trlIx ).*1e3;% convert spike times in s to ms
                                x(x<dt(1)) =[];
                                x(x>dt(end)) =[];
                                spkRaster(ft,:) = hist(x,dt);
                            end;
                            if size(spkRaster,1) ~= size(LogDat.LogDat1.log,1);
                                error('matrix dimension 1 must have equal size');
                            end;
                            
                            %% compute firing rate
                            T = bw;
                            step = bw*Fs;
                            FR{ct} = zeros(length(trl_all),3);
                            
                            tInt = [-1000:step:0];% work in ms
                            x = [];
                            for ft = 1:length(tInt)-1
                                nSpk = sum(spkRaster(:,find(dt >=tInt(ft) & dt < tInt(ft+1))),2);
                                x(:,ft) = nSpk./T;% firing rate
                            end;
                            FR{ct}(:,1) = mean(x,2);% mean firing rate during baseline
                            clear x;
                            
                            tInt = [500:step:1500];
                            x = [];
                            for ft = 1:length(tInt)-1
                                nSpk = sum(spkRaster(:,find(dt >=tInt(ft) & dt < tInt(ft+1))),2);
                                x(:,ft) = nSpk./T;% firing rate
                            end;
                            FR{ct}(:,2) = mean(x,2);% mean firing rate during cue
                            clear x;
                            
                            tInt = [2500:step:5000];
                            x = [];
                            for ft = 1:length(tInt)-1
                                nSpk = sum(spkRaster(:,find(dt >=tInt(ft) & dt < tInt(ft+1))),2);
                                x(:,ft) = nSpk./T;% firing rate
                            end;
                            FR{ct}(:,3) = mean(x,2);% mean firing rate during encoding
                            clear x;
                            
                            %% compute STATS
                            [~,~,~,stat1] = ttest(FR{ct}(:,2),FR{ct}(:,1));% compare spike count cue vs base
                            empT1(ct) = stat1.tstat;% save empirical t-stat
                            
                            [~,~,~,stat2] = ttest(FR{ct}(:,3),FR{ct}(:,1));% compare spike count encoding vs base
                            empT2(ct) = stat2.tstat;% save empirical t-stat
                            
                            design = [ones(1,nTrl(1)) 2*ones(1,nTrl(2)) 3*ones(1,nTrl(3))];
                            [empF1(ct),~,~] = computeFtest( design , FR{ct}(trlIx1,3));
                            
                            %design = [ones(1,nTrl(4)) 2*ones(1,nTrl(5)) 3*ones(1,nTrl(6))];
                            %[empF2(ct),~,~] = computeFtest( design , FR{ct}(trlIx2,3));
                            
                        end;
                    end;
                    fprintf('\n');
                end;
                
                %% create null-distribution from permutation data
                randT1 = [];        randT2 = [];
                randF1 = [];        randF2 = [];
                
                parfor pt = 1:nperms
                    %for pt = 1:nperms
                    
                    fprintf([num2str(pt),'/',num2str(nperms)]);
                    nTrl;
                    
                    rT1 = [];        rT2 = [];
                    rF1 = [];        rF2 = [];
                    
                    for it = 1:length(FR)%loop over all clusters
                        
                        x = [];    x = FR{it}(:,[1 2]);% base and cue
                        x2 = [];
                        for kt = 1:size(x,1)%loop over trials
                            x2(kt,:) = x(kt,randperm(size(x,2)));  % randomly swap the condition label
                        end;
                        
                        [~,~,~,stat1] = ttest(x2(:,2),x2(:,1));% compare spike count cue vs base
                        rT1(it) = stat1.tstat;% save empirical t-stat
                        
                        x =[];    x = FR{it}(:,[1 3]);% base and encoding
                        x2 = [];
                        
                        for kt = 1:size(x,1)%loop over trials
                            x2(kt,:) = x(kt,randperm(size(x,2))); % randomly swap the condition label
                        end;
                        
                        [~,~,~,stat2] = ttest(x2(:,2),x2(:,1));% compare spike count encoding vs base
                        rT2(it) = stat2.tstat;% save empirical t-stat
                        
                        design = [ones(1,nTrl(1)) 2*ones(1,nTrl(2)) 3*ones(1,nTrl(3))];
                        [rF1(it),~,~] = computeFtest( design , FR{ct}(randperm(size(FR{ct},1)),3));
                        
                        %design = [ones(1,nTrl(4)) 2*ones(1,nTrl(5)) 3*ones(1,nTrl(6))];
                        %[rF2(it),~,~] = computeFtest( design , FR{ct}(randperm(size(FR{ct},1)),3));
                        
                    end;% end of loop across clusters
                    
                    % keep the stats from this specific random partition
                    randT1(pt,:) = [min(rT1) max(rT1)];
                    randT2(pt,:) = [min(rT2) max(rT2)];
                    randF1(pt,:) = max(rF1);
                    %randF2(pt,:) = max(rF2);
                    
                    fprintf('\n');
                    
                end;% end of loop across permutations
                
                %% compute p-values
                pval1 = zeros(length(empT1),1);
                for it = 1:length(empT1);
                    if isnan(empT1(it));
                        pval1(it) = 1;
                    else
                        if sign(empT1(it))==1
                            pval1(it) = length(find(randT1(:,2)>=empT1(it)))/nperms;
                        else
                            pval1(it) = length(find(randT1(:,1)<=empT1(it)))/nperms;
                        end;
                    end;
                end;
                
                pval2 = zeros(length(empT2),1);
                for it = 1:length(empT2);
                    if isnan(empT2(it));
                        pval2(it) = 1;
                    else
                        if sign(empT2(it))==1
                            pval2(it) = length(find(randT2(:,2)>=empT2(it)))/nperms;
                        else
                            pval2(it) = length(find(randT2(:,1)<=empT2(it)))/nperms;
                        end;
                    end;
                end;
                
                pval3 = zeros(length(empF1),1);
                for it = 1:length(empF1);
                    if isnan(empF1(it));
                        pval3(it) = 1;
                    else
                        pval3(it) = length(find(randF1>=empF1(it)))/nperms;
                    end;
                end;
                
                %         pval4 = zeros(length(empF2),1);
                %         for it = 1:length(empF2);
                %             if isnan(empF2(it));
                %                 pval4(it) = 1;
                %             else
                %                 pval4(it) = length(find(randF2>=empF2(it)))/nperms;
                %             end;
                %         end;
                
                %%
                thr1 = zeros(size(clusterInfo,1),1);
                thr1(intersect(find(pval1 < 0.025),find(sign(empT1)==1)))  = 1;
                thr1(intersect(find(pval1 < 0.025),find(sign(empT1)==-1))) = -1;
                
                thr2 = zeros(size(clusterInfo,1),1);
                thr2(intersect(find(pval2 < 0.025),find(sign(empT2)==1)))  = 1;
                thr2(intersect(find(pval2 < 0.025),find(sign(empT2)==-1))) = -1;
                
                thr3 = zeros(size(clusterInfo,1),1);
                thr3(find(pval3 < 0.05)) = 1;
                
                %         thr4 = zeros(size(clusterInfo,1),1);
                %         thr4(find(pval4 < 0.05)) = 1;
                
                %% save the significance of effects for each cluster
                clusterInfo(:,11) = pval1;
                clusterInfo(:,13) = pval2;
                clusterInfo(:,15) = pval3;
                %clusterInfo(:,17) = pval4;
                
                clusterInfo(:,12) = thr1;
                clusterInfo(:,14) = thr2;
                clusterInfo(:,16) = thr3;
                %clusterInfo(:,18) = thr4;
                
                %%
                out_name = [savepath,'Cluster_info_EM_',pID{st},'_',sesh2 ,'.mat'];
                readme = {'sessionID','microwireID','clusterNR1','clusterNR2','initialClusterCategory','correctedClusterCategory','pctBelw','SNR','maxITC','peakSign','pval1','thr1','pval2','thr2','pval3','thr3'};
                %readme = {'sessionID','microwireID','clusterNR1','clusterNR2','initialClusterCategory','correctedClusterCategory','pctBelw','SNR','maxITC','peakSign','pval1','thr1','pval2','thr2','pval3','thr3','pval4','thr4'};
                save(out_name,'clusterInfo','CSClabel','CLUlabel','readme');
            catch
            end;
        end;
    end;
end;

%%
delete(gcp);
exit;





