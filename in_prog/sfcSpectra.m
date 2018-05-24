%%
restoredefaultpath;
addpath(genpath('~/tbx/chronux_2_11/'));

%%
pID = 'P04';
expMode = 'fVSpEM';
sesh = {'2016-10-23_15-08-47'};
seshIt = 1;

[p2LfpRes] = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',pID,'/res/'];
[p2SpkRes] = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',pID,'/spkRes/'];

[spectralFiles] = dir([p2LfpRes,pID,'_',expMode,'_',sesh{seshIt},'_*_spectralData.mat'])

%%
dt = -1e3:5e3;
load([p2LfpRes,pID,'_',expMode,'_',sesh{seshIt},'_preprocLFP.mat'])
[spkDat] = load([p2SpkRes,pID,'_',expMode,'_',sesh{seshIt},'_spkResp.mat'])

%%
T = 2;
W = 3;
TW = T*W;
k = 2*TW-1;

params1                  = [];
params1.Fs               = 1e3;
params1.fpass            = [0 30];
params1.pad              = 4;
params1.trialave         = 1;
params1.tapers           = [TW k];

res = 1/(2^(nextpow2(T*params1.Fs)+params1.pad)/params1.Fs);
[nfft1] = length( params1.fpass(1):res:params1.fpass(2) );

T = 2;
W = 10;
TW = T*W;
k = 2*TW-1;

params2                  = [];
params2.Fs               = 1e3;
params2.fpass            = [30 100];
params2.pad              = 4;
params2.trialave         = 1;
params2.tapers           = [TW k];

res = 1/(2^(nextpow2(T*params2.Fs)+params2.pad)/params2.Fs);
[nfft2] = length( params2.fpass(1)+res:res:params2.fpass(2) );

%%
[sfC1] = cell(length(spectralFiles),1);
[sfC2] = cell(length(spectralFiles),1);
[Sp1]  = cell(length(spectralFiles),1);
[Sp2]  = cell(length(spectralFiles),1);
[mSfc1] = zeros(length(spectralFiles),1);
[mSfc2] = zeros(length(spectralFiles),1);
[TFR1] = cell(length(spectralFiles),1);
[TFR2] = cell(length(spectralFiles),1);

%%
for curMW = 1:length(spectralFiles)
    
    fprintf([num2str(curMW),'/',num2str(length(lfpDat.LFPseg))]);
    [ntrl] = length(selTrl{curMW});
    
    [spectralDat] = load([p2LfpRes,spectralFiles(curMW).name]);
    startIx = regexp(spectralFiles(curMW).name,'\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}')+20;
    stopIx = regexp(spectralFiles(curMW).name,'_spectralData.mat')-1;
    [curMWname] = spectralFiles(curMW).name(startIx:stopIx);
    
    [~,pfIx] = max(spectralDat.S1);
    [pf] = spectralDat.f1(pfIx);
    
    M = ones(length(spectralDat.t3),1)*mean(spectralDat.S3,1);
    SD = ones(length(spectralDat.t3),1)*std(spectralDat.S3,0,1);
    TFR1{curMW} = (spectralDat.S3 - min(min(spectralDat.S3)))./(max(max(spectralDat.S3))-min(min(spectralDat.S3)));
    
    M = ones(length(spectralDat.t4),1)*mean(spectralDat.S4,1);
    SD = ones(length(spectralDat.t4),1)*std(spectralDat.S4,0,1);
    TFR2{curMW} = (spectralDat.S4 - min(min(spectralDat.S4)))./(max(max(spectralDat.S4))-min(min(spectralDat.S4)));
    
    %if (pf >= 5) && (pf <=12)
    
    % lfp data
    x1 = lfpDat.LFPseg{strcmp(lfpDat.chanLab,curMWname)};
    if isfield(lfpDat,'dsTrlTime')
        x1 = x1(lfpDat.dsTrlTime >= 2 & lfpDat.dsTrlTime <= 4,:);
    else
        x1 = x1(lfpDat.trlTime >= 2 & lfpDat.trlTime <= 4,:);
    end;
    
    % spike time data
    selIx = find(sum([spkDat.pval1 < 0.001 spkDat.pval2 < 0.001],2));
    
    %[sfC1{curMW}] = zeros(length(selIx),nfft1);
    %[sfC2{curMW}] = zeros(length(selIx),nfft2);
    %[Sp1{curMW}] = zeros(length(selIx),nfft1);    
    %[Sp2{curMW}] = zeros(length(selIx),nfft2);
    for curUnit = 1:length( selIx )
        x2 = spkDat.spkRaster{selIx(curUnit)}(:,dt>=2e3 & dt <=4e3)';
        ix = find(ismember(1:length(lfpDat.trlENC),selTrl{strcmp(lfpDat.chanLab,curMWname)}));
        x2 = x2(:,ix);
        
        if any(size(x1) ~= size(x2))
            error('lfp and spk matrix must have equal dimensions');
        end;
        
        %[sfC1{curMW}(curUnit,:),~,~,~,Sp1{curMW}(curUnit,:),f1] = coherencycpb(x1,x2,params1);
        %[sfC2{curMW}(curUnit,:),~,~,~,Sp2{curMW}(curUnit,:),f2] = coherencycpb(x1,x2,params2);
    end;
    %mSfc1(curMW) = max(max(sfC1{curMW}));
    %mSfc2(curMW) = max(max(sfC2{curMW}));
    %end;
    fprintf('\n');
end;

mx = [];
for it = 1:length(TFR1)
    mx(it) = max(max(TFR1{it}));
end;
outL = find(abs(mx-mean(mx))./std(mx)>3);
TFR1(outL) = [];

mx = [];
for it = 1:length(TFR2)
    mx(it) = max(max(TFR2{it}));
end;
outL = find(abs(mx-mean(mx))./std(mx)>3);
TFR2(outL) = [];

avgTFR = zeros(size(TFR1{1}));
for it = 1:length(TFR1)
    avgTFR = avgTFR+TFR1{it};
end;
avgTFR = avgTFR./it;
figure;
imagesc(spectralDat.t3,spectralDat.f3,avgTFR');
axis xy;

avgTFR = zeros(size(TFR2{1}));
for it = 1:length(TFR2)
    avgTFR = avgTFR+TFR2{it};
end;
avgTFR = avgTFR./it;
figure;
imagesc(spectralDat.t4,spectralDat.f4,avgTFR');
axis xy;

%%
figure;
for it = 1:length( stAVG)
    
    subplot(4,1,1:3);
    imagesc(stAVG{it});
    subplot(4,1,4);
    plot(-lag:lag,mean(stAVG{it},1));
    pause;
    clf;
    
end;

%%
cnt = 0;
        for curTrl = 1:ntrl
            [spkIx] = find( x2(:,curTrl) );            
            for curSpk = 1:length( spkIx )
                if (spkIx(curSpk)-lag>0) && (spkIx(curSpk)+lag<=size(x1,1))
                    cnt = cnt+1;
                    stAVG{curUnit}(cnt,:) = x1(spkIx(curSpk)-lag:spkIx(curSpk)+lag,curTrl);
                end;
            end;
        end;