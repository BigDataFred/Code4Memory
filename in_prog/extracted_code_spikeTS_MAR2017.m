    
%% load the Log data
    path_2_log = [rpath,sesh{gt},filesep,'log_dat',filesep];
    fn = dir([path_2_log,'*_EMtask_LogDat.mat']);
    
    LogDat = load([path_2_log,fn.name]);%load the data
    
    %% extract the trial information for the encoding phase
    [trl_all] = extract_Enc_trl(LogDat.LogDat1,LogDat.ix);% all, remembered (r), not-remembered (nr)
    

    %% load the LFP data for the selected seshion
    path_2_lfps = [rpath,sesh{gt},filesep,'lfp_dat',filesep];
    lfp_files = dir([path_2_lfps,'*.mat']);
    
    [lfp_dat] = load_lfp_data(lfp_files,path_2_lfps,trl_all);          
                        %% make the PSTH
                    n{ct} = zeros(length(dt),length(ts_trl));
                    for jt = 1:length(ts_trl)
                        n{ct}(:,jt) = hist(ts_trl{jt},dt);% count the number of spikes per time bin
                    end;
                    n{ct}(n{ct}>1)=1;
                    
                    sel_ix = find(dt >=-1e3 & dt <0);
                    [thr2] = mean(sum(n{ct}(sel_ix,:),1))+2*std(sum(n{ct}(sel_ix,:),1),0,2);
                    
                    %                     if floor(thr*1)/1 ~= floor(thr2*1)/1
                    %                         error('threshold estimates must be equal');
                    %                     end;
                    
                    %% smooth firing rate for visualization
                    nw = 50;
                    w = gausswin(21);%gaussian kernel
                    dt3 = -2e3:nw:4e3;% time bins with larger width
                    n2{ct} = zeros(length(dt3),length(ts_trl));
                    for jt = 1:length(ts_trl)
                        dum = hist(ts_trl{jt},dt3);
                        dum(dum>1)=1;
                        sn = length(dum);
                        dum = conv([fliplr(dum) dum fliplr(dum)],w,'same')./sum(w);%smooth the spike count over time
                        n2{ct}(:,jt) = dum(sn+1:2*sn);
                    end;
                    
                    %%
                    if ~isempty(n)% if the PSTH is not empty
                        
                        %% autocorrelation for each task interval
                        xcn = 1e3;
                        xc{ct} = zeros(xcn*2+1,size(n{ct},2));
                        xc1{ct} = zeros(xcn*2+1,size(n{ct},2));
                        xc2{ct} = zeros(xcn*2+1,size(n{ct},2));
                        xc3{ct} = zeros(xcn*2+1,size(n{ct},2));
                        for ft = 1:size(n{ct},2)
                            x = n{ct}(:,ft);%baseline
                            [xc{ct}(:,ft),lag] = xcorr(x-mean(x),xcn,'coeff');
                            x = n{ct}(i_base,ft);%baseline
                            [xc1{ct}(:,ft),lag] = xcorr(x-mean(x),xcn,'coeff');
                            x = n{ct}(i_cue,ft);%cue
                            [xc2{ct}(:,ft),lag] = xcorr(x-mean(x),xcn,'coeff');
                            x = n{ct}(i_enc,ft);%encoding
                            [xc3{ct}(:,ft),lag] = xcorr(x-mean(x),xcn,'coeff');
                        end;
                        
                        %% LFP data
                        cfg                     = [];
                        cfg.channel             = spike_dat{zt}.hdr.label;
                        
                        [dum] = ft_selectdata( cfg , lfp_dat );%select the relevant CSC channel
                        
                        % % put data into a conveniant format
                        cfg                     = [];
                        cfg.keeptrials          = 'yes';
                        
                        [dum] = ft_timelockanalysis( cfg, dum);
                        
                        %%
                        %perform interpolation of LFP in time domain to
                        %eliminate spike shape
                        int_step = [2 2];
                        for jt = 1:size(dum.trial,1)
                            lfp = squeeze( dum.trial(jt,:,:) );
                            spk = n{ct}(:,jt);
                            [lfp_interp] = interpLFP(lfp,spk,int_step);
                            
                            dum.trial(jt,:,:) = lfp_interp;
                        end;
                        
                        %extract LFP for different task intervals
                        lfp = squeeze(dum.trial)';
                        
                        lfp_base = lfp(i_base,:);
                        lfp_cue = lfp(i_cue,:);
                        lfp_enc = lfp(i_enc,:);
                        
                        
                    end;
%% STA
                        win = 250;
                        [STA{ct}] = compute_STA(win,n{ct},lfp);
                        
                        [STA1{ct}] = compute_STA(win,n{ct}(i_base,:),lfp_base);
                        [STA2{ct}] = compute_STA(win,n{ct}(i_cue,:),lfp_cue);
                        [STA3{ct}] = compute_STA(win,n{ct}(i_enc,:),lfp_enc);
                        
                        %% spectrum and SFC
                        TW = 3;
                        k = 2*TW-1;
                        
                        params                  = [];
                        params.Fs               = 1e3;
                        params.pad              = -1;
                        params.fpass            = [1 100];
                        params.tapers           = [TW k];
                        params.err              = [2 0.05];
                        params.trialave         = 1;
                        
                        [Snn1{ct},f1,R1,Serr1{ct}]=mtspectrumpb(n{ct}(i_base,:),params,0);
                        [Snn2{ct},f1,R2,Serr2{ct}]=mtspectrumpb(n{ct}(i_cue,:),params,0);
                        [Snn3{ct},f1,R3,Serr3{ct}]=mtspectrumpb(n{ct}(i_enc,:),params,0);
                        
                        %[~,dy] = gradient(lfp_base);
                        [Syy1{ct},f1,Serr4{ct}]=mtspectrumc(lfp_base,params);
                        %[~,dy] = gradient(lfp_cue);
                        [Syy2{ct},f1,Serr5{ct}]=mtspectrumc(lfp_cue,params);
                        %[~,dy] = gradient(lfp_enc);
                        [Syy3{ct},f1,Serr6{ct}]=mtspectrumc(lfp_enc,params);
                        
                        %                     [Syy1] = correct_one_over_f(Syy1, f1);
                        %                     [Syy2] = correct_one_over_f(Syy2, f1);
                        %                     [Syy3] = correct_one_over_f(Syy3, f1);
                        
                        TW = 12;
                        k = 2*TW-1;
                        
                        params                  = [];
                        params.Fs               = 1e3;
                        params.pad              = -1;
                        params.fpass            = [1 100];
                        params.tapers           = [TW k];
                        params.err              = [2 0.05];
                        params.trialave         = 1;
                        
                        [C1{ct},~,~,~,~,f2]=coherencycpb(lfp_base,n{ct}(i_base,:),params,0);
                        [C2{ct},~,~,~,~,f2]=coherencycpb(lfp_cue,n{ct}(i_cue,:),params,0);
                        [C3{ct},~,~,~,~,f2]=coherencycpb(lfp_enc,n{ct}(i_enc,:),params,0);
                        
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
                        [Cgrm{ct},~,~,~,~,tAx,fAx]=cohgramcpb(y,n{ct},movingwin,params,0);