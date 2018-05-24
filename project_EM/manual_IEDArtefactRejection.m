%% set libraries
addpath('/home/rouxf/tbx/fieldtrip-20170618/');
addpath('/home/rouxf/prj/Bham/code/mcode/helper/');
ft_defaults;

%% clear workspace
clear all;
clc;

%% set the main path
micPath = '/media/rouxf/rds-share/iEEG_DATA/MICRO/';
pDat = dir([micPath,'P*']);

%% loop over patients
for xt = 1:length(pDat)
    
    
    pID = pDat(xt).name;
    [rpath] = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',pID,'/'];
    savePath = rpath;
    
    %%
    expModes = dir(rpath);
    chck = regexp({expModes(:).name},'EM');
    sel = [];
    for it = 1:length(chck)
        
        if ~isempty(chck{it})
            if ( expModes(it).isdir )
                sel = [sel it];
            end;
        end;
    end;
    [expModes] = {expModes(sel).name};
    
    %% initialize
    dataAllLFPH = [];    dataAllLFPM = [];    dataAllSPKH = [];    dataAllSPKM = [];
    expLabH  = [];    expLabM  = [];    runLabH  = [];    runLabM  = [];
    cnt=0;    flag = 0;
    
    %% loop over expModes
    for yt = 1:length(expModes)
        
        [chck] = dir([savepath,pID,'preprocLFPdat_EMtask.mat']);
        
        if isempty(chck)
            rpath = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',pID,'/',expModes{yt},'/'];
            
            chck = dir(rpath);
            pattern = '\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}';
            sel = [];
            for it = 1:length(chck)
                if ~isempty(regexp(chck(it).name,pattern)) && (chck(it).isdir)
                    sel = [sel it];
                end
            end;
            sesh = {chck(sel).name}';
            
            %% loop over sessions
            for zt = 1:length(sesh)
                
                [datPath] = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',pID,'/',expModes{yt},'/',sesh{zt},'/'];
                
                %% search the Logfile
                p2logDat = [datPath,'log_dat/'];
                logFile = dir([p2logDat,'*_Log*ile_EMtask_LogDat.mat']);
                
                if ~isempty(logFile)
                    
                    %% load the logfile data
                    cnt = cnt+1;
                    
                    [logDat] = load([p2logDat,logFile.name]);
                    
                    [RTs]    = str2double(logDat.LogDat1.log(:,end));
                    [trlENC] = extract_Enc_trl(logDat.LogDat1);%1:length(logDat.LogDat1.log)+length(logDat.LogDat2.log);%
                    
                    %% indexes of trials that were later remebered and forgotten
                    [hitIdx] =logDat.ix{4};%HITS
                    [missIdx] = [logDat.ix{5};logDat.ix{6}];%MISSES
                    
                    %% delete trials with RTs < 4 sec
                    delIx = find(RTs(hitIdx) <=4);
                    hitIdx(delIx) = [];
                    
                    delIx = find(RTs(missIdx) <=4);
                    missIdx(delIx) = [];
                    
                    %% filter trial indexes
                    [trlENC1] = trlENC(hitIdx);
                    [trlENC2] = trlENC(missIdx);
                    
                    [RTs1]    = RTs(hitIdx);
                    [RTs2]    = RTs(missIdx);
                    
                    %% store exp and run information
                    [expLabH] = [expLabH yt*ones(1,length(trlENC1))];
                    [expLabM] = [expLabM yt*ones(1,length(trlENC2))];
                    
                    [runLabH] = [runLabH cnt*ones(1,length(trlENC1))];
                    [runLabM] = [runLabM cnt*ones(1,length(trlENC2))];
                    
                    %% read data
                    p2d = [datPath,'lfp_dat/'];
                    CSCfiles = dir([p2d,'*',sesh{zt},'_stimlocked.mat']);
                    
                    dum = cell(1,length(CSCfiles));
                    for it = 1:length(CSCfiles)
                        fprintf([num2str(it),'/',num2str(length(CSCfiles))]);
                        micDat = load([p2d,CSCfiles(it).name]);
                        dum{it} = micDat.save_data{1}{1}{1};
                        fprintf('\n');
                    end;
                    
                    [dumLFP] = ft_appenddata( [] , dum{:} );
                    dumLFP.fsample = dumLFP.cfg.previous{1}.resamplefs;
                    
                    %% 50 Hz cleaning
                    signal = dumLFP.trial;
                    dumLFP.trial{1} = {};
                    for it = 1:length(dumLFP.label)
                        for jt = 1:length(signal)
                            [cleanSignal,~] = CleanLineNoise( signal{jt}(it,:) ,'Fs', dumLFP.fsample , 'noiseFreq', 50,'windowSize',1);
                            signal{jt}(it,:) = cleanSignal;
                        end;
                        cleanSignal = [];
                    end;
                    dumLFP.trial = signal;
                    signal = [];
                    
                    %% split data into hits and misses
                    cfg                     = [];
                    cfg.trials              = trlENC1;
                    
                    [dumLFP1] = ft_selectdata( cfg , dumLFP );
                    
                    cfg                     = [];
                    cfg.trials              = trlENC2;
                    
                    [dumLFP2] = ft_selectdata( cfg , dumLFP );
                    clear dumLFP;
                    
                    %% read data
                    p2d = [datPath,'spike_dat/'];
                    SPKfiles = dir([p2d,'*',sesh{zt},'_stimlocked.mat']);
                    dumSPK1 = cell(1,length(SPKfiles));
                    dumSPK2 = cell(1,length(SPKfiles));
                    for it = 1:length(SPKfiles)
                        fprintf([num2str(it),'/',num2str(length(SPKfiles))]);
                        spkDat = load([p2d,SPKfiles(it).name]);
                        dumSPK = spkDat.save_data{1}{1}{1};
                        
                        for jt = 1:length(dumSPK.label)
                            dumSPK.label(jt) = {[dumSPK.hdr.label,'_',dumSPK.label{jt},num2str(cnt)]};
                        end;
                        
                        dumSPK1{it} = dumSPK;
                        dumSPK2{it} = dumSPK;
                        for jt = 1:length(dumSPK.unit)
                            if ~isempty(dumSPK.unit{jt})
                                [selIx1] = find(ismember(dumSPK.trial{jt},trlENC1));
                                [selIx2] = find(ismember(dumSPK.trial{jt},trlENC2));
                                
                                dumSPK1{it}.trial{jt} = dumSPK1{it}.trial{jt}(selIx1);
                                dumSPK1{it}.timestamp{jt} = dumSPK1{it}.timestamp{jt}(selIx1);
                                dumSPK1{it}.time{jt} = dumSPK1{it}.time{jt}(selIx1);
                                dumSPK1{it}.unit{jt} = dumSPK1{it}.unit{jt}(selIx1);
                                dumSPK1{it}.waveform{jt} = dumSPK1{it}.waveform{jt}(:,:,selIx1);
                                
                                dumSPK2{it}.trial{jt} = dumSPK2{it}.trial{jt}(selIx2);
                                dumSPK2{it}.timestamp{jt} = dumSPK2{it}.timestamp{jt}(selIx2);
                                dumSPK2{it}.time{jt} = dumSPK2{it}.time{jt}(selIx2);
                                dumSPK2{it}.unit{jt} = dumSPK2{it}.unit{jt}(selIx2);
                                dumSPK2{it}.waveform{jt} = dumSPK2{it}.waveform{jt}(:,:,selIx2);
                            else
                                dumSPK1{it}.trial{jt} = [];
                                dumSPK1{it}.timestamp{jt} = [];
                                dumSPK1{it}.time{jt} = [];
                                dumSPK1{it}.unit{jt} = [];
                                dumSPK1{it}.waveform{jt} = [];
                                
                                dumSPK2{it}.trial{jt} = [];
                                dumSPK2{it}.timestamp{jt} = [];
                                dumSPK2{it}.time{jt} = [];
                                dumSPK2{it}.unit{jt} = [];
                                dumSPK2{it}.waveform{jt} = [];
                            end;
                        end;
                        fprintf('\n');
                    end;
                    clear dumSPK;
                    
                    [dumSPK1] = ft_appendspike(  [] , dumSPK1{:});
                    [dumSPK2] = ft_appendspike(  [] , dumSPK2{:});
                    
                    %% concetenate LFP-data across recording sessions
                    if cnt == 1
                        if ~isempty(trlENC1) > 0
                            [dataAllLFPH] = dumLFP1;
                        end;
                        if ~isempty(trlENC2) > 0
                            [dataAllLFPM] = dumLFP2;
                        end;
                        [hRTs] = RTs1;
                        [mRTs] = RTs2;
                    else
                        if ~isempty(trlENC1) > 0
                            [dataAllLFPH] = ft_appenddata([], dataAllLFPH, dumLFP1);
                        end;
                        if ~isempty(trlENC2) > 0
                            [dataAllLFPM] = ft_appenddata([], dataAllLFPM, dumLFP2);
                        end;
                        [hRTs] = [hRTs;RTs1];
                        [mRTs] = [mRTs;RTs2];
                    end;
                    clear dumLFP1 dumLFP2;
                    
                    %% concetenate SPK-data across recording sessions
                    if cnt == 1
                        if ~isempty(trlENC1) > 0
                            [dataAllSPKH] = dumSPK1;
                        end;
                        if ~isempty(trlENC2) > 0
                            [dataAllSPKM] = dumSPK2;
                        end;
                        [hRTs] = RTs1;
                        [mRTs] = RTs2;
                    else
                        if ~isempty(trlENC1) > 0
                            [dataAllSPKH] = ft_appendspike([], dataAllSPKH, dumSPK1);
                        end;
                        if ~isempty(trlENC2) > 0
                            [dataAllSPKM] = ft_appendspike([], dataAllSPKM, dumSPK2);
                        end;
                        [hRTs] = [hRTs;RTs1];
                        [mRTs] = [mRTs;RTs2];
                    end;
                    clear dumSPK1 dumSPK2;
                    
                end;
                
            end;
        else
            load([savepath,chck.name]);
            flag = 1;
        end;
    end;
    
    %% save data concatenated over sessions
    if flag ==0
        saveName = [pID,'_concatenatedLFPdat_EMtask.mat'];
        save([savePath,saveName],'dataAllLFPH','dataAllLFPM','hRTs','mRTs','expLabH','expLabM','runLabH','runLabM');
        
        saveName = [pID,'_concatenatedSPKdat_EMtask.mat'];
        save([savePath,saveName],'dataAllSPKH','dataAllSPKM','hRTs','mRTs','expLabH','expLabM','runLabH','runLabM');
    end;
            
    %% mark IED artefacts by hand
    cfg                     = [];
    cfg.lpfilter            = 'yes';
    cfg.lpfreq              = 100;
    cfg.demean              = 'yes';
    cfg.detrend             = 'yes';
    
    [dum1] = ft_preprocessing( cfg , dataAllLFPH );
    
    cfg                     = [];
    cfg.latency             = [-.5 4.5];
    
    [dum1] = ft_selectdata( cfg , dum1 );
    
    for it = 1:length( dum1.trial )
        fprintf([num2str(it),'/',num2str(length( dum1.trial ))]);
        x = dum1.trial{it};
        M  = mean(x,2);
        SD = std(x,0,2);
        M  = M*ones(1,size(dum1.trial{it},2));
        SD = SD*ones(1,size(dum1.trial{it},2));
    
        dum1.trial{it} = (x-M)./SD;
        fprintf('\n');
    end;
    clear x M SD;
    
    cfg                     = [];
    cfg.viewmode            = 'vertical';
    
    [cfg1] = ft_databrowser( cfg, dum1);
    
    %% mark IED artefacts by hand
    cfg                     = [];
    cfg.lpfilter            = 'yes';
    cfg.lpfreq              = 100;
    cfg.demean              = 'yes';
    cfg.detrend             = 'yes';
    
    [dum2] = ft_preprocessing( cfg , dataAllLFPM );
    
    cfg                     = [];
    cfg.latency             = [-.5 4.5];
    
    [dum2] = ft_selectdata( cfg , dum2 );
    
    for it = 1:length( dum2.trial )
        fprintf([num2str(it),'/',num2str(length( dum2.trial ))]);
        x = dum2.trial{it};
        M  = mean(x,2);
        SD = std(x,0,2);
        M  = M*ones(1,size(dum2.trial{it},2));
        SD = SD*ones(1,size(dum2.trial{it},2));
    
        dum2.trial{it} = (x-M)./SD;
        fprintf('\n');
    end;
    clear x M SD;
    
    cfg                     = [];
    cfg.viewmode            = 'vertical';
    
    [cfg2] = ft_databrowser( cfg, dum2);
    
    %% return to command line
    return;
    clc;
    
    %% set bad channels manually
    badChans = {'CSC_antHippR1','CSC_antHippR2','CSC_midHippR7','CSC_midHippR8','CSC_postHippR4','CSC_postHippR5','CSC_postHippR6','CSC_postHippR8'};
    [goodChans] = setdiff( dataAllLFPH.label , badChans );
    
    %% create dummy smaples
    trl = [];
    ix = 1:length( dum1.trial{1});
    for it = 1:length(dum1.trial)
    
        trl(it,:) = ix;
        ix = ix + length(dum1.trial{it});
    
    end;
    clear ix;
    
    %% extract trials with artefacts
    delIx = cfg1.artfctdef.visual.artifact;
    chck = [];
    for jt = 1:size(trl,1)
        x = [];
        for it = 1:length(delIx)
    
            if any(ismember(trl(jt,:),delIx(it,:)))
                x = [x it];
            end;
        end;
        chck{jt} = x;
    
    end;
    
    [BADtrlIdx1] = [];
    c=0;
    for it = 1:length( chck )
    
        if ~isempty(chck{it})
            c = c+1;
            BADtrlIdx1= [BADtrlIdx1 it];
        end;
    
    end;
    
    %% create dummy smaples
    trl = [];
    ix = 1:length( dum2.trial{1});
    for it = 1:length(dum2.trial)
    
        trl(it,:) = ix;
        ix = ix + length(dum2.trial{it});
    
    end;
    clear ix;
    
    %% extract trials with artefacts
    delIx = cfg2.artfctdef.visual.artifact;
    chck = [];
    for jt = 1:size(trl,1)
        x = [];
        for it = 1:length(delIx)
    
            if any(ismember(trl(jt,:),delIx(it,:)))
                x = [x it];
            end;
        end;
        chck{jt} = x;
    
    end;
    
    [BADtrlIdx2] = [];
    c=0;
    for it = 1:length( chck )
    
        if ~isempty(chck{it})
            c = c+1;
            BADtrlIdx2= [BADtrlIdx2 it];
        end;
    
    end;
    
    %% extract artefact-free trials
    [GOODtrlIdx1] = setdiff(1:length(dataAllLFPH.trial),BADtrlIdx1);
    [GOODtrlIdx2] = setdiff(1:length(dataAllLFPM.trial),BADtrlIdx2);
    
    %% save the good and bad trial indexes
    saveName = [pID,'_IEDdetection_trialIndexes.mat'];
    save([savePath,saveName],'BADtrlIdx1','BADtrlIdx2','GOODtrlIdx1','GOODtrlIdx2','goodChans','badChans');
    
    %% trim RTs and label data
    hRTs(BADtrlIdx1) = [];
    mRTs(BADtrlIdx2) = [];
    expLabH(BADtrlIdx1) = [];
    expLabM(BADtrlIdx2) = [];
    runLabH(BADtrlIdx1) = [];
    runLabM(BADtrlIdx2) = [];
    
    %% trim the SPK data
    cfg                     = [];
    cfg.trials              = GOODtrlIdx1;
    cfg.channel             = goodChans;
    
    [dataAllLFPH] = ft_selectdata( cfg , dataAllLFPH );
    
    cfg                     = [];
    cfg.trials              = GOODtrlIdx2;
    cfg.channel             = goodChans;
    
    [dataAllLFPM] = ft_selectdata( cfg , dataAllLFPM );
    
    %% trim the SPK data 
    for jt = 1:length(dataAllSPKH.unit)
        if ~isempty(dataAllSPKH.unit{jt})                        
            [selIx1] = find(ismember(dataAllSPKH.trial{jt},GOODtrlIdx1));
            [selIx2] = find(ismember(dataAllSPKM.trial{jt},GOODtrlIdx2));
            
            dataAllSPKH.trial{jt} = dataAllSPKH.trial{jt}(selIx1);
            dataAllSPKH.timestamp{jt} = dataAllSPKH.timestamp{jt}(selIx1);
            dataAllSPKH.time{jt} = dataAllSPKH.time{jt}(selIx1);
            dataAllSPKH.unit{jt} = dataAllSPKH.unit{jt}(selIx1);
            dataAllSPKH.waveform{jt} = dataAllSPKH.waveform{jt}(:,:,selIx1);
            
            dataAllSPKM.trial{jt} = dataAllSPKM.trial{jt}(selIx2);
            dataAllSPKM.timestamp{jt} = dataAllSPKM.timestamp{jt}(selIx2);
            dataAllSPKM.time{jt} = dataAllSPKM.time{jt}(selIx2);
            dataAllSPKM.unit{jt} = dataAllSPKM.unit{jt}(selIx2);
            dataAllSPKM.waveform{jt} = dataAllSPKM.waveform{jt}(:,:,selIx2);
        else
            dataAllSPKH.trial{jt} = [];
            dataAllSPKH.timestamp{jt} = [];
            dataAllSPKH.time{jt} = [];
            dataAllSPKH.unit{jt} = [];
            dataAllSPKH.waveform{jt} = [];
            
            dataAllSPKM.trial{jt} = [];
            dataAllSPKM.timestamp{jt} = [];
            dataAllSPKM.time{jt} = [];
            dataAllSPKM.unit{jt} = [];
            dataAllSPKM.waveform{jt} = [];
        end;
    end;
     
    %%   
    delIx = [];
    for jt = 1:length(dataAllSPKH.unit)
        
        [chanLab] = dataAllSPKH.label{jt}(1:regexp(dataAllSPKH.label{jt},'_sig')-1);
        
        if (isempty(dataAllSPKH.trial{jt})) || ( ismember(chanLab,badChans) ) || ( length(dataAllSPKH.timestamp{jt}) < 50 ) || ( unique(dataAllSPKH.unit{jt}) ==0 )
            delIx = [delIx jt];
        end;
    end;
    
    dataAllSPKH.label(delIx) = [];
    dataAllSPKH.trial(delIx) = [];
    dataAllSPKH.timestamp(delIx) = [];
    dataAllSPKH.time(delIx) = [];
    dataAllSPKH.unit(delIx) = [];
    dataAllSPKH.waveform(delIx) = [];
    
    dataAllSPKM.label(delIx) = [];
    dataAllSPKM.trial(delIx) = [];
    dataAllSPKM.timestamp(delIx) = [];
    dataAllSPKM.time(delIx) = [];
    dataAllSPKM.unit(delIx) = [];
    dataAllSPKM.waveform(delIx) = [];
                
    %% save the trimmed data
    saveName = [pID,'_preprocLFPdat_EMtask.mat'];
    save([savePath,saveName],'dataAllLFPH','dataAllLFPM','hRTs','mRTs','expLabH','expLabM','runLabH','runLabM','-v7.3');
    saveName = [pID,'_preprocSPKdat_EMtask.mat'];
    save([savePath,saveName],'dataAllSPKH','dataAllSPKM','hRTs','mRTs','expLabH','expLabM','runLabH','runLabM','-v7.3');
    
end;
