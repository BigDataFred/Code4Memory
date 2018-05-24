function computeSpkLFPparameters4EMtask(dataID,nCPUs)

tic;

if nargin == 0
    [  dataID  ] = 'P07'; % IDs can be P02, P04 , P05, P07, P08, P22AMS, P23AMS
    nCPUs = 36;
end;

%% set the path environment
addpath('/home/rouxf/tbx/fieldtrip-20170618/');
ft_defaults;
addpath(genpath('/home/rouxf/tbx/chronux_2_11/'));
addpath('/home/rouxf/tbx/chronux_2_11/spectral_analysis/statistical_tests/');
addpath('/home/rouxf/prj/Bham/code/mcode/utils/LFP/');
addpath('/home/rouxf/prj/Bham/code/mcode/utils/SFC/');
addpath('/home/rouxf/prj/Bham/code/mcode/project_EM/');
addpath('/home/rouxf/prj/Bham/code/mcode/utils/spikes/');
addpath('/home/rouxf/prj/Bham/code/mcode/utils/hybrid/');
addpath('/home/rouxf/prj/Bham/code/mcode/helper/');
addpath('/home/rouxf/prj/Bham/code/mcode/visualize/continuous/');
addpath('/home/rouxf/tbx/CircStat2012a/');

%%
if isempty(gcp('NoCreate'))
    parpool(nCPUs,'SpmdEnabled',false);
end

%% parameters controlling which data to load
[expMode] = {'fvSpEM' 'cnEM'};% experiment modes

[MWlabel] = {'1' '2' '3' '4' '5' '6' '7' '8' '*'};% labels of the MWs
bfSel      = [];% default = [] (all), to select BFs  use 1-nBFs
bfSel2     = bfSel;% default = [] (all),to select BFs use 1-nBFs
mwSel      = 9; % select MW (1,2,3,4,5,6,7,8,9), note that 9 is equal to 'all'
micDatMode = 0;% SFC analysis with single MW (0), locally re-referenced MW (1), average MW on BF (2)

coreVars = {'dataID','expMode','MWlabel','bfSel','bfSel2','mwSel','micDatMode','nCPUs','coreVars','patIt','exIt','seshIt','expSel',...
    'sesh','BFlabel','seshLab','dataMode','BFlabel2'};% never erase these variables

%%
for patIt = 1:length(dataID)% loop over patients
    
    for exIt = 1:length(expMode) % loop over experiments
        
        expSel     = exIt;% switch between fVSpEM (1) & cnEM (2)
        
        %%
        if iscell(dataID)
            pID = dataID{patIt};
        else
            pID = dataID;
        end;
        
        %% extract session and BG labels
        [sesh,BFlabel,seshLab,dataMode] = extractSessionandBFlabels(pID,expMode,expSel,[]);
        BFlabel2 = BFlabel;
        
        for seshIt = 1:length(sesh)% loop over sessions
            
            %%
            if iscell(dataID)
                pID = dataID{patIt};
            else
                pID = dataID;
            end;
            
            %%
            try
                seshSel    = seshIt;% switch date and time of recording (1,nsesh)
                
                %% set savepath and savename for results
                savepath = '/media/rouxf/rds-share/iEEG_DATA/MICRO/code4MEM/';
                savename = [pID,'_',expMode{expSel},'_',sesh{seshSel}];
                %savename = [pID,'_',expMode{expSel},'_',sesh{seshSel},'_pooledDataSPKandLFP.mat'];
                
                %% set up path and filenames of Micro, Macro and LogData
                switch dataMode
                    case 'BHX'
                        [p2SpkDat,spkFiles,p2micDat,micFiles,p2macDat,macFile,macDat,ixIEDdat,logDat] = setupMicANDMacANDLogPathANDFileNames4EMtask(pID,expMode,expSel,sesh,seshSel,BFlabel,BFlabel2,bfSel,bfSel2,MWlabel,mwSel,micDatMode);
                    case 'AMS'
                        [p2SpkDat,spkFiles,p2micDat,micFiles,p2macDat,macFile,macDat,ixIEDdat,logDat] = setupMicANDMacANDLogPathANDFileNames4EMtaskAMS(pID,expMode,expSel,sesh,seshSel,BFlabel,BFlabel2,bfSel,bfSel2,MWlabel,mwSel,micDatMode);
                end;
                %%
                savelist = {'hitIdx','missIdx','hitIdx2','missIdx2','trlENC','BFinf','glFR','FR','spkRaster','uSelIx','trXC','instFR','wvfStats','XC','data_lfp','data_lfp2',...
                    'data_spk', 'data_spk2','LFPr','phi','pvalPL','hPL','medLFPsig','spectrumPOW','allSTP1','STApow1','allSTP2','STApow2','spectrumSFC','spectralData','staDat','LFP2LFPcoh','chID','ITC','PPC1','PPC2','savepath','savename'};
                
                %% extract the trials corresponding to the encoding
                [trlENC] = extract_Enc_trl(logDat.LogDat1);%1:length(logDat.LogDat1.log)+length(logDat.LogDat2.log);%
                
                %% indexes of trials that were later remebered and forgotten
                [hitIdx] =logDat.ix{4};%HITS
                [missIdx] = [logDat.ix{5};logDat.ix{6}];%MISSES
                trlENC(missIdx) = [];
                
                %% read the spk data
                [spkDat,uIx,spkRaster,chLab] = readSpkdat(p2SpkDat,spkFiles,trlENC);
                
                %% filter the clusters and visualize Units
                plt = 'yes';
                %[spkDat,BFinf,glFR,FR,uSelIx,trXC,instFR,wvfStats] = clusterFilt(x,spkDat,uIx,trlENC,seshLab,plt);
                [spkDat,BFinf,glFR,FR,uSelIx,instFR,wvfStats,spkRaster] = clusterFilt(spkRaster,spkDat,uIx,trlENC,seshLab,plt);
                trlENCspk = trlENC;
                
                save([savepath,savename,'_spikeANDwvfDataAll.mat'],'BFinf','glFR','FR','uSelIx','instFR','wvfStats','spkRaster','trlENC','-v7.3');
                                
                %% compute the auto and crosscorrelation matrix for selected units
                plt = 'no';
                [XC] = computeClusterXcorr(spkDat,[2 5].*1e3,480,trlENC,plt,chLab);
                
                %%
                selIdx = 3;
                trXC = cell(1,length(selIdx));
                tt = [-1 5].*1e3;
                tt = tt(1):tt(2);
                for kt = 1:length(selIdx)
                    [trXC{kt}.xc,trXC{kt}.fr,trXC{kt}.lag,trXC{kt}.ntrl] = time_resolved_xcorr(tt,480,trlENC,spkDat{selIdx(kt)},[]);
                end;
                    
                %% import the MW-LFP-data
                [micDat] = importMICdat(p2micDat,micFiles);
                
                %% delete those trial labels that have been flaged by IED detection
                [trlENC,IEDidx,MicTrlIx1,MicTrlIx2] = deleteIEDflaggedTRL(micDat,ixIEDdat,macDat,trlENC);
                
                %                 %% indexes of trials that were later remebered and forgotten
                %                 [hitIdx2] = find(ismember(trlENC,hitIdx));%HITS
                %                 [missIdx2] = find(ismember(trlENC,missIdx));%MISSES
                %                 if length(hitIdx2)+length(missIdx2) ~= length(trlENC)
                %                     error('trial indexes must be consistent');
                %                 end;
                
                %% extract the trials that survived LFP-artefact rejection from spike data
                [spkDat] = trlSelectspkDat(spkDat,trlENC);
                
                %% extract the MW-LFP data corresponding to the Encoding trials
                [origMWdat,IEDdat,chID,origTime] = extractMWandIEDdat(micDat,micDatMode,MicTrlIx1,MicTrlIx2,IEDidx);
                [dt] = origMWdat.time{1}(2)-origMWdat.time{1}(1);
                Fs = 1/dt;
                
                %FIXME
                [BFinf] = recodeBFinf(BFinf,origMWdat.label);
                
                %% search for sneaky IEDs
                [trl,pct,chLabIED] = IEDdetector(origMWdat,BFinf(1,:),[-1 5]);
                delChan = find(pct >=.35);
                delIxC = {};
                for it = 1:length(delChan)                    
                    delIxC{it} = find(strcmp(BFinf(1,:),chLabIED(delChan(it))));
                end;
                delIxC = [delIxC{:}];
                
                keepChan = setdiff(1:length(BFinf(1,:)),delIxC);
                
                delIx = trl(setdiff(1:length(trl),delChan));
                delIx = unique([delIx{:}]);                
                trlENC(delIx) = [];                                
                
                if ~isempty(trlENC) && length(trlENC) >= 20 && ~isempty(keepChan)
                    
                    cfg             = [];
                    cfg.trials      = setdiff(1:length(origMWdat.trial),delIx);
                    cfg.channels    = BFinf(1,keepChan);
                    
                    [origMWdat] = ft_selectdata( cfg ,origMWdat );
                    
                    spkDat = spkDat(keepChan);
                    BFinf = BFinf(:,keepChan);
                    
                    %% do local referencing of MW-LFP- data
                    if micDatMode ==1 % locally re-reference the MWs to the weakest MW
                        [origMWdat] = localRefMWdat(origMWdat);
                    end;
                    
                    %% load the Macro-LFP-data and extract trials corresponding to Encoding
                    if ~isempty( macDat )
                        [MACdat] = loadMACdat(micDat,macDat,trlENC);
                    end;
                    
                    %% demean and detrend
                    [origMWdat] = preprocLFP(origMWdat);
                    
                    cfg                     = [];
                    cfg.bpfilter            = 'yes';
                    cfg.bpfreq              = [3 8];
                    cfg.bpfilttype          = 'but';
                    cfg.bpfiltord           = 4;
                    
                    [BPdat] = preprocLFP(origMWdat,cfg);
                    
                    %% converts data from ft to chx format (inludes interp of LFP at spike times)
                    [data_lfp ,data_spk, data_spk2] = ft2chronux(spkDat,origMWdat,BFinf,trlENC,[min(origMWdat.time{1}) max(origMWdat.time{1})]);
                    [data_lfpBP] = ft2chronux([],BPdat,BFinf,trlENC,[min(origMWdat.time{1}) max(origMWdat.time{1})]);
                    chID = chID([BFinf{2,:}]);
                    
                    %% make a sanity check of the LFP label assignments
                    sanityCheckLFP2SPKassigment(origMWdat,spkDat,BFinf,data_lfp);
                    
                    %% check number of spikes across conditions, this makes sure that all units have at least equal minimum number of spikes per cond
                    minSpk = 50;
                    
                    selIdx = find(origTime >= -1 & origTime <=5);
                    
                    [chck] = checkSpikeNumberPerCond([],[],selIdx,data_spk,minSpk,1);% each unit must fire at least x*2 spikes
                    
                    %% only keep data of units that satisfy minimum number of spikes
                    selIdx = find(chck);
                    inList = {'spkDat','BFinf','glFR','FR','instFR', 'wvfStats', 'spkRaster', 'data_lfp','data_lfpBP','data_spk','data_spk2','chID'};
                    [outDat] = selectChannelData(selIdx,{spkDat,BFinf,glFR, FR,instFR, wvfStats, spkRaster,data_lfp, data_lfpBP, data_spk,data_spk2,chID});
                    
                    uSelIx = uSelIx(selIdx,:);
                    spkRaster = spkRaster(selIdx,:);
                    
                    if length(inList) == length(outDat)
                        for it = 1:length(inList)
                            eval([inList{it},'= outDat{',num2str(it),'};']);
                        end;
                    else
                        error('number of output and input arguments must have equal length');
                    end;
                    
                    dat4Phi = data_lfp;
                    
                    save([savepath,savename,'_spikeANDwvfData.mat'],'BFinf','glFR','FR','uSelIx','instFR','wvfStats','spkRaster','trlENC','trlENCspk','-v7.3');
                    
                    %%
                    [ITC,ERPimg,tAx] = computeITC4LFP(origMWdat,[1:30],[-1 5],BFinf(1,:));
                    
                    save([savepath,savename,'_ITC.mat'],'ITC','ERPimg','tAx','BFinf');
                    
                    %                 %% compute the pairwise correlation matrix
                    %                 id = unique({BFinf{1,:}});
                    %                 ixMED = zeros(1,length(id));
                    %                 for it = 1:length(id)
                    %                     ixMED(it) = min(find(strcmp(BFinf(1,:),id(it))));
                    %                 end;
                    %
                    %                 plt = 'n';
                    %                 [LFPr] = computeBFcorrLFP_datModeA(data_lfp(ixMED),'Spearman',plt);%FIXME: also implement LFP-LFP coherence
                    %                 x = sign(LFPr(:));
                    %                 if any(x ==-1)
                    %                     %     error('Phase jump detected');
                    %                 end;
                    
                    %% extract phase
                    foi = 1:30;
                    toi = origTime;
                    %[phi,phiIx] = extractPhaseLFP(MWdat,foi,toi);
                    selIdx = find(origTime >= toi(1) & origTime <=toi(end));
                    dum = cell(1,length(dat4Phi));
                    for it = 1:length(dat4Phi)
                        dum{it} = dat4Phi{it}(selIdx,:);
                    end;
                    [phi,phiIx] = extractPhaseLFP(dum,foi,toi,Fs);
                    
                    selIdx = find(ismember(phiIx,find(origTime >= 2 & origTime <=4)));
                    for it = 1:length( phi )
                        phi{it} = phi{it}(:,:,selIdx);
                    end;
                    clear dum;
                    
%                     %% extract time range of interest for lfp and spike data
%                     selIdx = find(origMWdat.time{1} >= -1 & origMWdat.time{1} <=5);
%                     [data_lfp] = extractTimeRangeOfInterestChronuxData(selIdx,[-1 5],data_lfp);
%                     [data_lfpBP] = extractTimeRangeOfInterestChronuxData(selIdx,[-1 5],data_lfpBP);
%                     [data_spk] = extractTimeRangeOfInterestChronuxData(selIdx,[-1 5],data_spk);
%                     [data_spk2] = extractTimeRangeOfInterestChronuxData(selIdx,[-1 5],data_spk2);
%                     
%                     cfg                     = [];
%                     cfg.latency             = [-1 5] ;
%                     
%                     [MWdat] = ft_selectdata(cfg , origMWdat);
%                     
%                     %% do pre-whitening (1/f correction)
%                     pwMethod = 'dx/dt';%'dx/dt'; % can be dx/dt or ARMfilt
%                     [data_lfp2] = prewhitenLFP(data_lfp,pwMethod,Fs);
%                     
%                     save([savepath,savename,'_LFPandSpkData.mat'],'data_lfp','data_spk','data_spk2','data_lfpBP','data_lfp2','BFinf','-v7.3');
%                     
%                     %%
%                     timeGrid = MWdat.time{1};
%                     
%                     PPC = zeros(length(data_lfp),30);
%                     for it = 1:length(data_lfp)
%                         [PPC(it,:),freqAx] = computePPC(data_lfp{it},data_spk{it},Fs,timeGrid,'n');
%                     end;
%                     save([savepath,savename,'_PPC.mat'],'PPC','freqAx','-v7.3');
                    
                    %% compute spike-phase locking
                    foi = 1:30;
                    pvalPL = cell(1,length(data_spk));
                    hPL = zeros(1,length(data_spk));
                    zPL= cell(1,length(data_spk));
                    PL= cell(1,length(data_spk));
                    for it = 1:length(data_spk)
                        %[pvalPL{it},zPL{it}] = computeSpkPL(data_spk{it}(phiIx,:),phi{BFinf{2,it}});
                        [pvalPL{it},zPL{it},PL{it}] = computeSpkPL(data_spk{it},phi{it});
                        hPL(it) = any(pvalPL{it} < (0.05/ length(foi)));
                    end;
                    save([savepath,savename,'_PLV.mat'],'hPL','pvalPL','zPL','PL','foi','phi','-v7.3');
                    
%                     %%
%                     [spectrumPOW,spectrumSFC] = computeSpectrumAndSFC4timeRange(data_lfp , data_spk, 1:size(data_lfp{1},1),1:size(data_lfp{1},2), Fs);
%                     save([savepath,savename,'_spectrumPOWandSFC.mat'],'spectrumPOW','spectrumSFC','-v7.3');
%                     
%                     %%
%                     smp =MWdat.time{1};
%                     D = [-.5 .5];
%                     plt = 'n';
%                     err = 0;
%                     allSTA = cell(1,length(data_lfp));
%                     allSTAbp = cell(1,length(data_lfp));
%                     allSTP = cell(1,length(data_lfp));
%                     STApow = cell(1,length(data_lfp));
%                     for it = 1:length(data_lfp)
%                         [allSTA{it}]    = sta(data_spk2{it},data_lfp{it},smp,plt,[],[-1 5],D,err);
%                         [allSTAbp{it}]    = sta(data_spk2{it},data_lfpBP{it},smp,plt,[],[-1 5],D,err);
%                         [allSTP{it}]    = computeSpectrumAndSFC4timeRange({allSTA{it}'}, [], [], [],Fs);
%                         [STApow{it}]    = computeSpectrumAndSFC4timeRange({mean(allSTA{it},1)'}, [], [], [],Fs);
%                     end;
%                     save([savepath,savename,'_STAandSTPandSTApow.mat'],'allSTA','allSTAbp','allSTP','STApow','-v7.3');
%                     
%                     smp =MWdat.time{1};
%                     D = [-.5 .5];
%                     plt = 'n';
%                     err = 0;
%                     allSTP1 = cell(1,length(data_lfp));
%                     STApow1 = cell(1,length(data_lfp));
%                     allSTP2 = cell(1,length(data_lfp));
%                     STApow2 = cell(1,length(data_lfp));
%                     for it = 1:length(data_lfp)
%                         [allSTA]    = sta(data_spk2{it},data_lfp{it},smp,plt,[],[0 2],D,err);
%                         if size(allSTA,2) >1
%                             [allSTP1{it}]    = computeSpectrumAndSFC4timeRange({allSTA'}, [], [], [],Fs);
%                             [STApow1{it}]    = computeSpectrumAndSFC4timeRange({mean(allSTA,1)'}, [], [], [],Fs);
%                         end;
%                         [allSTA]    = sta(data_spk2{it},data_lfp{it},smp,plt,[],[2 4],D,err);
%                         if size(allSTA,2)>1
%                             [allSTP2{it}]    = computeSpectrumAndSFC4timeRange({allSTA'}, [], [], [],Fs);
%                             [STApow2{it}]    = computeSpectrumAndSFC4timeRange({mean(allSTA,1)'}, [], [], [],Fs);
%                         end;
%                     end;
%                     save([savepath,savename,'_STPandSTApowCueANDEncoding.mat'],'allSTP1','STApow1','allSTP2','STApow2','-v7.3');
                    
                    %                 %%
                    %                 if length(pvalPL) ~= length(data_lfp);error('data dimensions out of range');end;
                    %
                    %                 staDat = cell(1,sum(hPL));
                    %                 cnt = 0;
                    %                 for kt = 1:length(data_lfp)
                    %                     fprintf([num2str(kt),'/',num2str(length(data_lfp))])
                    %                     if hPL(kt) ~=0
                    %                         cnt = cnt+1;
                    %                         [staDat{cnt}] = computeSFC4HitsAndMisses(data_lfp{kt},data_lfpBP{kt},data_spk2{kt},hitIdx2,missIdx2,MWdat.time{1},Fs);
                    %                     end;
                    %                     fprintf('\n');
                    %                 end;
                    
%                     %% Time-Frequency analysis
%                     [spectralData] = runSpectralAnalysis4LFPandSPK(data_lfp, data_lfp, data_spk,Fs);
%                     save([savepath,savename,'_timeFrequencyData.mat'],'spectralData','-v7.3');
%                     
                    %                 %% compute the auto and crosscorrelation matrix for selected units
                    %                 plt = 'no';
                    %                 [XC] = computeClusterXcorr(spkDat(hPL==1),[-1 5].*1e3,480,trlENC,plt,chLab);
                    %
                    %                 selIdx = find(hPL==1);
                    %                 trXC = cell(1,length(selIdx));
                    %                 tt = [-1 5].*1e3;
                    %                 tt = tt(1):tt(2);
                    %                 for kt = 1:length(selIdx)
                    %                     [trXC{kt}.xc,trXC{kt}.fr,trXC{kt}.lag,trXC{kt}.ntrl] = time_resolved_xcorr(tt,480,trlENC,spkDat{selIdx(kt)},[]);
                    %                 end;
                    %
                    %                 %% extract the median LFP - use for phase estimate
                    %                 [medLFPsig] = extractMedianLFP(data_lfp(ixMED),LFPr,chID(ixMED));% median of raw signal
                    %                 %[medLFPsig2] = extractMedianLFP(data_lfp2(ixMED),LFPr,chID(ixMED));% median of pre-whitened signal
                    %
                    %                 %%
                    %                 [LFP2FLPcoh] = computeLFP2LFPcoh(medLFPsig,Fs);
                    %
                    %                 %%
                    %                 [LFPr] = squeeze(mean(LFPr,1));
                    
                    %                 %% clear useless data and make estimate of disc size
                    %                 x = whos;
                    %                 wsp = {};
                    %                 for it = 1:length(x)
                    %                     wsp  = [wsp x(it).name];
                    %                 end;
                    %                 clearList = {setdiff(wsp,savelist)};
                    %                 clearList = {setdiff(clearList{1},coreVars),'wsp','x','clearList'};
                    %                 delIdx = [];
                    %                 for it = 1:length(clearList{1})
                    %                     if any(strcmp(clearList{1}(it),coreVars))
                    %                         cnt = cnt+1;
                    %                         delIdx(cnt) = it;
                    %                     end;
                    %                 end;
                    %                 if isempty(delIdx)
                    %                     clearList = [clearList{:}];
                    %                     for it = 1:length(clearList);clear(clearList{it});end;
                    %                 end;
                    %                 x = whos;
                    %                 B = [];
                    %                 for it = 1:length(x)
                    %                     B(it) = x(it).bytes;
                    %                 end;
                    %
                    %                 diskSpace = sum(B)/1e9;
                    %                 clear B;
                    
                    %                 %% save data to disc
                    %                 if diskSpace > 2
                    %                     save([savepath,savename],'-v7.3');
                    %                 else
                    %                     save([savepath,savename]);
                    %                 end;
                    toc;
                end;
            catch ME
                savename = [pID,'_',expMode{expSel},'_',sesh{seshSel},'_pooledDataSPKandLFP.mat'];
                [~,savename,~] = fileparts(savename);
                errFile = ['errorLog_',savename,'.txt'];
                fid = fopen([savepath,errFile],'w');
                errorMSG = ME.message;
                fprintf(fid,errorMSG);
                fprintf(fid,'\n');
                fclose(fid);
            end;
        end;% end of loop over sessions
    end;% end of loop over experimetns
end;% end of loop over patients


return;