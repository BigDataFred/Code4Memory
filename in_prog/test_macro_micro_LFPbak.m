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
    parpool(36,'SpmdEnabled',false);
end
tic;

%% parameters controlling which data to load
[  pID  ] = 'P07'; % IDs can be P02, P04 , P05, P07, P08, P22AMS, P23AMS
[expMode] = {'fvSpEM' 'cnEM'};% experiment modes
[MWlabel] = {'1' '2' '3' '4' '5' '6' '7' '8' '*'};% labels of the MWs

expSel     = 1;% switch between fVSpEM (1) & cnEM (2)
seshSel    = 1;% switch date and time of recording (1,nsesh)
bfSel      = [];% default = [] (all), to select BFs  use 1-nBFs
bfSel2     = bfSel;% default = [] (all),to select BFs use 1-nBFs
mwSel      = 9; % select MW (1,2,3,4,5,6,7,8,9), note that 9 is equal to 'all'
micDatMode = 0;% SFC analysis with single MW (0), locally re-referenced MW (1), average MW on BF (2)

%% extract session and BG labels
[sesh,BFlabel,seshLab,dataMode] = extractSessionandBFlabels(pID,expMode,expSel,seshSel);
BFlabel2 = BFlabel;

%% set savepath and savename for results
savepath = '/media/rouxf/rds-share/iEEG_DATA/MICRO/code4MEM/';
savename = [pID,'_',expMode{expSel},'_',sesh{seshSel},'_pooledDataSPKandLFP.mat'];

%% set up path and filenames of Micro, Macro and LogData
switch dataMode
    case 'BHX'
        [p2SpkDat,spkFiles,p2micDat,micFiles,p2macDat,macFile,macDat,ixIEDdat,logDat] = setupMicANDMacANDLogPathANDFileNames4EMtask(pID,expMode,expSel,sesh,seshSel,BFlabel,BFlabel2,bfSel,bfSel2,MWlabel,mwSel,micDatMode);
    case 'AMS'
        [p2SpkDat,spkFiles,p2micDat,micFiles,p2macDat,macFile,macDat,ixIEDdat,logDat] = setupMicANDMacANDLogPathANDFileNames4EMtaskAMS(pID,expMode,expSel,sesh,seshSel,BFlabel,BFlabel2,bfSel,bfSel2,MWlabel,mwSel,micDatMode);
end;
%%
savelist = {'hitIdx','missIdx','trlENC','BFinf','glFR','FR','uSelIx','trXC','instFR','wvfStats','XC','data_lfp','data_lfp2',...
    'data_spk', 'data_spk2','LFPr','phi','pvalPL','hPL','medLFPsig','spectrumPOW','spectrumSFC','hitsSFC','missSFC','hitsSTA',...
    'missSTA','spectralData','hitsSTP','missSTP','hitsSTAbp','missSTAbp','savepath','savename'};

%% extract the trials corresponding to the encoding
[trlENC] = extract_Enc_trl(logDat.LogDat1);%1:length(logDat.LogDat1.log)+length(logDat.LogDat2.log);%

%% indexes of trials that were later remebered and forgotten
[hitIdx] = trlENC(logDat.ix{4});%HITS
[missIdx] = trlENC([logDat.ix{5};logDat.ix{6}]);%MISSES

%% read the spk data
[spkDat,uIx,x,chLab] = readSpkdat(p2SpkDat,spkFiles,trlENC);

%% filter the clusters and visualize Units
plt = 'no';
%[spkDat,BFinf,glFR,FR,uSelIx,trXC,instFR,wvfStats] = clusterFilt(x,spkDat,uIx,trlENC,seshLab,plt);
[spkDat,BFinf,glFR,FR,uSelIx,instFR,wvfStats] = clusterFilt(x,spkDat,uIx,trlENC,seshLab,plt);

%% import the MW-LFP-data
[micDat] = importMICdat(p2micDat,micFiles);

%% delete those trial labels that have been flaged by IED detection
[trlENC,IEDidx,MicTrlIx1,MicTrlIx2] = deleteIEDflaggedTRL(micDat,ixIEDdat,macDat,trlENC);

%% indexes of trials that were later remebered and forgotten
[hitIdx] = find(ismember(trlENC,hitIdx));%HITS
[missIdx] = find(ismember(trlENC,missIdx));%MISSES

%% extract the trials that survived LFP-artefact rejection from spike data
[spkDat] = trlSelectspkDat(spkDat,trlENC);

%% extract the MW-LFP data corresponding to the Encoding trials
[origMWdat,IEDdat,chID,origTime] = extractMWandIEDdat(micDat,micDatMode,MicTrlIx1,MicTrlIx2,IEDidx);
[dt] = origMWdat.time{1}(2)-origMWdat.time{1}(1);
Fs = 1/dt;

%FIXME
[BFinf] = recodeBFinf(BFinf,origMWdat.label);

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
cfg.bpfreq              = [2 8];
cfg.bpfilttype          = 'but';
cfg.bpfiltord           = 4;

[BPdat] = preprocLFP(origMWdat);

%% converts data from ft to chx format (inludes interp of LFP at spike times)
[data_lfp ,data_spk, data_spk2] = ft2chronux(spkDat,origMWdat,BFinf,trlENC,[min(origMWdat.time{1}) max(origMWdat.time{1})]);
[data_lfpBP] = ft2chronux([],BPdat,BFinf,trlENC,[min(origMWdat.time{1}) max(origMWdat.time{1})]);
chID = chID([BFinf{2,:}]);

%% make a sanity check of the LFP label assignments
sanityCheckLFP2SPKassigment(origMWdat,spkDat,BFinf,data_lfp);

%% check number of spikes across conditions, this makes sure that all units have at least equal minimum number of spikes per cond
minSpk = 50;
[chck] = checkSpikeNumberPerCond(hitIdx,missIdx,data_spk,minSpk);

%% only keep data of units that satisfy minimum number of spikes in both conditions
selIdx = find(chck);
inList = {'spkDat','BFinf','glFR','FR','uSelIx','instFR', 'wvfStats', 'data_lfp','data_lfpBP','data_spk','data_spk2','chID'};
[outDat] = selectChannelData(selIdx,{spkDat,BFinf,glFR, FR,uSelIx,instFR, wvfStats, data_lfp, data_lfpBP, data_spk,data_spk2,chID});

if length(inList) == length(outDat)
    for it = 1:length(inList)
        eval([inList{it},'= outDat{',num2str(it),'};']);
    end;        
else
    error('number of output and input arguments must have equal length');
end;

dat4Phi = data_lfp;

%% compute the pairwise correlation matrix
id = unique({BFinf{1,:}});
ixMED = zeros(1,length(id));
for it = 1:length(id)
    ixMED(it) = min(find(strcmp(BFinf(1,:),id(it))));
end;

plt = 'n';
[LFPr] = computeBFcorrLFP_datModeA(data_lfp(ixMED),'Spearman',plt);%FIXME: also implement LFP-LFP coherence
x = sign(LFPr(:));
if any(x ==-1)
%     error('Phase jump detected');
end;

%% extract phase
foi = 2.^([6:2:64]./8);
toi = origTime;
%[phi,phiIx] = extractPhaseLFP(MWdat,foi,toi);
selIdx = find(origTime >= toi(1) & origTime <=toi(end));
dum = cell(1,length(dat4Phi));
for it = 1:length(dat4Phi)
    dum{it} = dat4Phi{it}(selIdx,:);
end;
[phi,phiIx] = extractPhaseLFP(dum,foi,toi,Fs);

selIdx = find(ismember(phiIx,find(origTime >= -1 & origTime <=5)));
for it = 1:length( phi )
    phi{it} = phi{it}(:,:,selIdx);
end;
clear dum;

%% extract time range of interest for lfp and spike data
selIdx = find(origMWdat.time{1} >= -1 & origMWdat.time{1} <=5);
[data_lfp] = extractTimeRangeOfInterestChronuxData(selIdx,[-1 5],data_lfp);
[data_lfpBP] = extractTimeRangeOfInterestChronuxData(selIdx,[-1 5],data_lfpBP);
[data_spk] = extractTimeRangeOfInterestChronuxData(selIdx,[-1 5],data_spk);
[data_spk2] = extractTimeRangeOfInterestChronuxData(selIdx,[-1 5],data_spk2);

cfg                     = [];
cfg.latency             = [-1 5] ;

[MWdat] = ft_selectdata(cfg , origMWdat);

%% compute spike-phase locking
pvalPL = cell(1,length(data_spk));
hPL = zeros(1,length(data_spk));
zPL= cell(1,length(data_spk));
for it = 1:length(data_spk)
    %[pvalPL{it},zPL{it}] = computeSpkPL(data_spk{it}(phiIx,:),phi{BFinf{2,it}});
    [pvalPL{it},zPL{it}] = computeSpkPL(data_spk{it},phi{it});
    hPL(it) = any(pvalPL{it}(foi >=2 & foi <=8) < 0.05/6);
end;

%% do pre-whitening (1/f correction)
pwMethod = 'dx/dt';%'dx/dt'; % can be dx/dt or ARMfilt
[data_lfp2] = prewhitenLFP(data_lfp,pwMethod,Fs);

%%
[spectrumPOW,spectrumSFC] = computeSpectrumAndSFC4timeRange(data_lfp , data_spk, 1:size(data_lfp{1},1),1:size(data_lfp{1},2), Fs);

%%
if length(pvalPL) ~= length(data_lfp);error('data dimensions out of range');end;

hitsSFC = cell(1,sum(hPL));
missSFC = cell(1,sum(hPL));
hitsSTA = cell(1,sum(hPL));
missSTA = cell(1,sum(hPL));
hitsSTP = cell(1,sum(hPL));
missSTP = cell(1,sum(hPL));
hitsSTAbp = cell(1,sum(hPL));
missSTAbp = cell(1,sum(hPL));
cnt = 0;
for kt = 1:length(data_lfp)
    fprintf([num2str(kt),'/',num2str(length(data_lfp))])
    if hPL(kt) ~=0
        cnt = cnt+1;
        [hitsSFC{cnt},missSFC{cnt},hitsSTA{cnt},missSTA{cnt},hitsSTP{cnt},missSTP{cnt},hitsSTAbp{cnt},missSTAbp{cnt}] = computeSFC4HitsAndMisses(data_lfp{kt},data_lfpBP{kt},data_spk2{kt},hitIdx,missIdx,MWdat.time{1},Fs);        
    end;
    fprintf('\n');
end;

%% Time-Frequency analysis
[spectralData] = runSpectralAnalysis4LFPandSPK(data_lfp(hPL==1), data_lfp2(hPL==1), data_spk(hPL==1),Fs);

%% compute the auto and crosscorrelation matrix for selected units
plt = 'no';
[XC] = computeClusterXcorr(spkDat(hPL==1),[-1 5].*1e3,480,trlENC,plt,chLab);

selIdx = find(hPL==1);
trXC = cell(1,length(selIdx));
tt = [-1 5].*1e3;
tt = tt(1):tt(2);
for kt = 1:length(selIdx)
    [trXC{kt}.xc,trXC{kt}.fr,trXC{kt}.lag,trXC{kt}.ntrl] = time_resolved_xcorr(tt,480,trlENC,spkDat{selIdx(kt)},[]);
end;

%% extract the median LFP - use for phase estimate
[medLFPsig] = extractMedianLFP(data_lfp(ixMED),LFPr,chID(ixMED));% median of raw signal
%[medLFPsig2] = extractMedianLFP(data_lfp2(ixMED),LFPr,chID(ixMED));% median of pre-whitened signal

LFPr = squeeze(mean(LFPr,1));

%% clear useless data and make estimate of disc size
x = whos;
wsp = {};
for it = 1:length(x)
    wsp  = [wsp x(it).name];
end;
clearList = {setdiff(wsp,savelist),'savelist'};
clearList = [clearList{:}];
for it = 1:length(clearList);clear(clearList{it});end;
clear clearList wsp x;

x = whos;
B = [];
for it = 1:length(x)
    B(it) = x(it).bytes;
end;

diskSpace = sum(B)/1e9;
clear B;

%% save data to disc
if diskSpace > 2
    save([savepath,savename],'-v7.3');
else
    save([savepath,savename]);
end;
toc;

return;