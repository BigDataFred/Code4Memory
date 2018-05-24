%%
restoredefaultpath;
addpath('~rouxf/prj/Bham/code/mcode/project_EM/');
addpath(genpath('~rouxf/prj/Bham/code/mcode/utils/'));
addpath(genpath('~rouxf/prj/Bham/code/mcode/visualize/'));
addpath(genpath('~/tbx/chronux_2_11/spectral_analysis/'));
addpath('~/tbx/eeglab14_1_1b/functions/sigprocfunc/');
addpath('~/tbx/fieldtrip-20170618/');
ft_defaults;

%% recruit workers for parallel computing
if isempty(gcp('nocreate'))
    parpool(36,'SpmdEnabled', false)
end;
tic;

%%
Fs = 1e3;
[hpF] = fir1(3*floor(Fs/0.5),0.5./(Fs/2),'high');

% t = 0:1/Fs:14;
% sig = sin(2*pi*0.5.*t)+0.5*sin(2*pi*2.*t)+.25*sin(2*pi*15.*t)+.25*randn(1,length(t));
% filt = filtfilt(hpF,1,[sig fliplr(sig) sig fliplr(sig) sig]);
% filt = filt(2*length(sig)+1:length(sig)*3);
% %filt2 = eegfilt(sig,Fs,0.5,0);
% figure;
% hold on;
% plot(t,sig);
% %plot(t,filt2,'r');
% plot(t,filt,'g');
% xlim([0 4]);

%%
[bpF] = fir1(3*fix(Fs/3),[3 8]./(Fs/2),'bandpass');

% t = 0:1/Fs:30;
% sig = .25*sin(2*pi*8.*t)+.25*randn(1,length(t));
% filt = filtfilt(bpF,1,[fliplr(sig) sig fliplr(sig)]);
% filt = filt(length(sig)+1:length(sig)*2);
% filt2 = eegfilt(sig,Fs,4,11);
% figure;
% hold on;
% plot(t,sig);
% plot(t,filt2,'r');
% plot(t,filt,'g');
% xlim([0 1]);

%%
[pID] = {'P02','P04','P05','P07','P08','P09'};%

[expMode] = {'fVSpEM','cnEM'};

%%
seshCnt = 0;
chanCntH = 0;
chanCntM = 0;
trlCntH  = 0;
trlCntM  = 0;
dt = -1500:5500;
dt2 = -1500:250:5500;
filtERP = [];
avgLFPh = [];
LFPpoolH = {};
avgLFPm = [];
LFPpoolM = {};
iFRpoolH = [];
iFRpoolM = [];
iFRpool = [];
unitPoolLabel = {};
lfpPoolLabel = {};
spkPool = {};
frPool = [];
delIxPool = {};
ixHPool = {};
ixMPool = {};

%%
for curPat = 1:length(pID)
    for curExp  = 1:length( expMode )
        %%
        
        %p2d = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',pID{curPat},'/',expMode{curExp},'/'];
        p2d = ['/home/rouxf/tmp/resultsC4M/',pID{curPat},'/'];
        
        %sesh = dir(p2d);
        sesh = dir([p2d,pID{curPat},'_',expMode{curExp},'*_lfpDataStimLockedSegmentedAVGdownsampled.mat']);
        
        if ~isempty(sesh)
            %sesh(1:2) = [];
            chck = regexp({sesh(:).name},'\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}');
            %selIx= []; cnt = 0;for it = 1:length(chck); if (~isempty(chck{it})) && (sesh(it).isdir);cnt = cnt+1; selIx(cnt) = it;end;end;
            selIx= []; cnt = 0;for it = 1:length(chck); if (~isempty(chck{it}));cnt = cnt+1; selIx(cnt) = it;end;end;
            %sesh = {sesh(selIx).name}'
            dum = cell(1,length(sesh));
            for it = 1:length(chck)
                dum(it) = {sesh(selIx(it)).name(chck{it}:chck{it}+18)}'
            end;
            sesh = dum
            
            %%
            for curSesh = 1:length( sesh )
                
                seshCnt = seshCnt +1;
                
                %%
                spkFile = [ pID{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_spkDataStimLockedSegmented.mat' ];
                %[spkDat] = load([p2d,sesh{curSesh},'/',spkFile])
                [spkDat] = load([p2d,spkFile])
                
                logFile = [ pID{curPat},'_',sesh{curSesh},'_LogFile_EMtask_LogDat.mat' ];
                %load([p2d,sesh{curSesh},'/',logFile])
                
                [trlIx] = spkDat.trlENC;
                [hitIdx] = find(ismember(sort([spkDat.hitIdx;spkDat.missIdx]),spkDat.hitIdx));
                [missIdx] = find(ismember(sort([spkDat.hitIdx;spkDat.missIdx]),spkDat.missIdx));
                
                try
                    load([p2d,logFile])
                    [stimCat] = logDat.LogDat1.log(:,3:4);
                    stimEv = cell(1,size(stimCat,1));
                    for stimIt = 1:size( stimCat ,1 )
                        stimEv(stimIt) = {[stimCat{stimIt,1}(1) stimCat{stimIt,2}(1)]};
                    end;
                    stimID = unique(stimEv);                                        
                catch     
                    [trlIx] = 1:length(spkDat.trlENC);
                    stimID = {};
                    stimID(1) = {'dum'};
                    stimEv = cell(1,length(trlIx));
                    for it = 1:length( trlIx )
                        stimEv(it) = stimID;
                    end;
                end;                                
                
                stimEv = stimEv(sort([hitIdx;missIdx]));
                    
                [cluDat] = extractClusterDat2( spkDat );
                
                %%
                sigSel = [];
                if length(hitIdx) >= 20
                    for curClust = 1:length( cluDat.chanLab )
                        
                        [spkRaster] = zeros(length(trlIx),length(dt));
                        [iFR] = zeros(length(trlIx),length(dt));
                        [fr] = zeros(length(trlIx),length(dt2));
                        parfor curTrl = 1:length(trlIx)
                            x = cluDat.ts{curClust}(cluDat.trl{curClust} == trlIx(curTrl));
                            [n,~] = hist(x,dt);
                            spkRaster(curTrl,:) = n;
                            iFR(curTrl,:) = conv(n,gausswin(201),'same')./sum(gausswin(201));
                            [n,~] = hist(x,dt2);
                            fr(curTrl,:) = n;
                        end;
                        [fr] = fr(:,dt2>=-1000 & dt2 <=5e3);
                        fr = mean(fr,1)./0.25;
                        tx2 = dt2(dt2>=-1000& dt2 <=5e3);
                        
                        [iFR] = iFR(:,dt>=-1000 & dt <=5e3);
                        [spkRaster] = spkRaster(:,dt>=-1000 & dt <=5e3);
                        tx = dt(dt>=-1000 & dt <=5e3);
                        
                        [iFR] = iFR./1e-3;
                        trsh = mean(mean(iFR(:,tx<=0),2))+1.5*std(mean(iFR(:,tx<=0),1));
                        if trsh < 1
                            trsh = 1;
                        end;
                         
                        iFRH = [];iFRM = [];
                        iFRH = mean(iFR(hitIdx,:),1);
                        if ~isempty(missIdx)
                            iFRM = mean(iFR(missIdx,:),1);
                        end;
                        iFR = mean(iFR,1);
                        
                        xix1 = find(sign(diff(iFR(tx>=200)>trsh))==1)+1;
                        xix2 = find(sign(diff(iFR(tx>=200)>trsh))==-1)+1;

                        if ~isempty(xix1) || ~isempty(xix2)
                            if min(xix2) < min(xix1)
                                xix1 = [1 xix1];
                            end;
                            if length(xix1)>length(xix2)
                                xix2(end+1) = length(iFRH);
                            end;
                        end;
                        
                        [baseCnt] = sum(spkRaster(:,tx>=-400 & tx <=0),2);
                        baseSD = 2*std(baseCnt);
                        
                        baseLine = median(baseCnt);
                        
                        if any(diff([xix1' xix2'],[],2)>100)                                                        

                            pval = [];
                            for curStim = 1:length(stimID)
                                stimIx = find(strcmp(stimEv,stimID{curStim}));
                                TW = [200 600];
                                p = zeros(1,44);
                                h = zeros(1,44);                               
                                for curTW = 1:44
                                    [postCnt] = sum(spkRaster(:,tx>=TW(1) & tx <=TW(2)),2);
                                    if median(postCnt(stimIx))>0
                                        [p(curTW)] = ranksum(postCnt(stimIx),baseCnt(stimIx));
                                        %[statInfo] = baseVSpostDepTtest4SPKdat(baseCnt,postCnt);
                                        TW = TW+100;
                                        h(curTW) = 1;
                                    end;
                                end;
                                pval = [pval min(p(h==1))];
                            end;
                            
                            if any(pval<0.01) 
                                sigSel = [sigSel curClust];
                                frPool = [frPool;fr];
                                iFRpool = [iFRpool;iFR];
                                iFRpoolH = [iFRpoolH;iFRH];
                                unitPoolLabel = [unitPoolLabel [pID{curPat},'-',sesh{curSesh},'-',cluDat.chanLab{curClust}]];
                                spkPool = [spkPool spkRaster];
                                
                                if ~isempty( iFRM)
                                    iFRpoolM = [iFRpoolM;iFRM];
                                end;
                                
                                %visualizeSpkRaster(tx,spkRaster,min(pval),tx2,iFR,trsh,fr); 

                            end;
                        end;
                        
                    end;
                end;
                
                %%
                if ~isempty( sigSel)
                    
                    sigSel = unique(sigSel);
                    
                    chanIx = [];
                    for curChan = 1:length(sigSel)
                        [chanIx(curChan)] = find(strcmp(spkDat.chanLab,cluDat.chanLab(sigSel(curChan))));
                    end;
                    
                    chanIx = unique(chanIx);                                                       
                    
                    %%
                    LFPfn = [pID{curPat},'_',expMode{curExp},'_',sesh{curSesh},'_lfpDataStimLockedSegmenteddownsampled.mat'];
                    [~,fn,ext] = fileparts(LFPfn);
                    
                    %[lfpDat] = load( [p2d,sesh{curSesh},'/',LFPfn] )
                    [lfpDat] = load( [p2d,LFPfn] )
                    if ( ~isfield(lfpDat,'dsTrlTime') )
                        lfpDat.dsTrlTime = lfpDat.trlTime;
                    end;
                    
                    %%
                    IEDfn = [fn,'_IEDdetection',ext];
                    %load([p2d,sesh{curSesh},'/',IEDfn]);
                    load([p2d,IEDfn]);
                    
                    %%
                    for jt = 1:length(chanIx)
                        
                        fprintf([num2str(jt),'/',num2str(length(chanIx))]);
                        
                        [lfp] = lfpDat.LFPseg{chanIx(jt)}(:,1:length(lfpDat.trlENC));
                        [trlIx] = 1:length(lfpDat.trlENC);
                        [trlPool] = sort([lfpDat.hitIdx;lfpDat.missIdx]);
                        
                        dum = spkPool{end-(length(chanIx)-(length(chanIx)-(length(chanIx)-jt)))};                        
                        if ismember(chanIx(jt),selIx{1})
                            lfp(:,[delIx{1}]) = [];
                            trlIx([delIx{1}]) = [];
                            trlPool([delIx{1}]) = [];
                            dum([delIx{1}],:) = [];
                        else
                            lfp(:,[delIx{2}]) = [];
                            trlIx([delIx{2}]) = [];
                            trlPool([delIx{2}]) = [];
                            dum([delIx{2}],:) = [];
                        end;
                        
                        [ixH] = find(ismember(trlPool,lfpDat.hitIdx));
                        [ixM] = find(ismember(trlPool,lfpDat.missIdx));
                        
                        dum(ixH,:);
                        
                        tIx = find(lfpDat.dsTrlTime >=-.5 & lfpDat.dsTrlTime <=-.1);
                        %M = ones(size(lfp))*mean(mean(lfp(tIx,:),1));
                        %SD = ones(size(lfp))*std(mean(lfp(tIx,:),1));
                        M = ones(size(lfp,1),1)*mean(lfp(tIx,:),1);
                        SD = ones(size(lfp,1),1)*std(lfp(tIx,:),1);
                        lfp = (lfp-M)./SD;
                        
                        tIx = find(lfpDat.dsTrlTime >=-1 & lfpDat.dsTrlTime <=5);
                        if length(ixH) >= 20
                            
                            if ismember(chanIx(jt),selIx{1})
                                if ~isempty(delIx{1})
                                    delIxPool = [delIxPool delIx{1}];
                                else
                                    delIxPool = [delIxPool 0];
                                end;
                            else
                                if ~isempty(delIx{2})
                                    delIxPool = [delIxPool delIx{2}];
                                else
                                    delIxPool = [delIxPool 0];
                                end;
                            end;
                            
                            ixHPool = [ixHPool ixH];                            
                            
                            lfpPoolLabel = [lfpPoolLabel [pID{curPat},'-',sesh{curSesh},'-',spkDat.chanLab{chanIx(jt)}]];
                                               
                            lfpH = lfp(:,ixH);
                            chanCntH = chanCntH+1;
                            trlCntH = trlCntH + length(ixH);                            
                            
                            avgLFPh = [avgLFPh mean(lfp(tIx,:),2)];                           
                            LFPpoolH = [LFPpoolH lfpH(tIx,:)];                            
                            
                            if ~isempty(ixM)
                                lfpM = lfp(:,ixM);
                                chanCntM = chanCntM+1;
                                trlCntM = trlCntM + length(ixM);                                                                
                                avgLFPm = [avgLFPm mean(lfpM(tIx,:),2)];                                
                                LFPpoolM = [LFPpoolM lfpM(tIx,:)];
                                ixMPool = [ixMPool ixM];
                            else
                                avgLFPm = [avgLFPm NaN(size(mean(lfpH(tIx,:),2)))];
                                LFPpoolM = [LFPpoolM NaN(size(lfpH(tIx,:)))];
                                ixMPool = [ixMPool 0];
                            end;
                        end;
                        
                        fprintf('\n');
                    end;
                    
                end;
            end;
        end;
    end;
end;
toc;

%% spike-interpolation step of LFP
for curLFP = 1:length(lfpPoolLabel)
    
    [lfp] = LFPpoolH{curLFP};
    [nsmp,ntrl] = size(lfp);
    
    unitIx = find(strcmp(lfpPoolLabel(curLFP),unitPoolLabel));
    
    for curUnit = 1:length(unitIx)
        spk = spkPool{unitIx(curUnit)};
        if (delIxPool{curLFP}~=0)
            spk(delIxPool{curLFP},:) = [];
        end;
        spk = spk(ixHPool{curLFP},:)';
        
        for curTrl = 1:ntrl
            tsIx = find(spk(:,curTrl)==1);
            if ~isempty(tsIx)
                lfp(:,curTrl) = interpLFP(lfp(:,curTrl),tsIx,1e3,'linear');
            end;
        end;
        LFPpoolH{curLFP} = lfp;
    end;
end;

%% 1/f correction of LFP
% pwLFP = cell(1,length(LFPpoolH));
% for it = 1:length( LFPpoolH)
%     pwLFP(it) = prewhitenLFP(LFPpoolH(it),'ARMfilt',1e3);
% end;

%%
iFRzsc = zeros(size(iFRpool));
for it = 1:size(iFRpool,1)
    x = iFRpool(it,:);
    M = nanmean(x(tx > -500 & tx <=0))*ones(1,size(x,2));
    SD = nanstd(x(tx > -500 & tx <=0))*ones(1,size(x,2));
    iFRzsc(it,:) = (x-M)./SD;
end;

%%
frZsc = zeros(size(frPool));
for it = 1:size(frPool,1)
    x = frPool(it,:);
    M = nanmean(x(tx2 > -500 & tx2 <=0))*ones(1,size(x,2));
    SD = nanstd(x(tx2 > -500 & tx2 <=0))*ones(1,size(x,2));
    frZsc(it,:) = (x-M)./SD;
end;
frZsc(find(isinf(frZsc(:,1))),:) = [];

uIx1 = [];
uIx2 = [];
for it = 1:size(frZsc,1)
    
    [~,ix] = max(frZsc(it,:));
    if tx2(ix) > 0 && tx2(ix) < 2e3
        uIx1 = [uIx1 it];
    elseif tx2(ix) > 2e3
        uIx2 = [uIx2 it];
    end;
    
end;

figure;
hold on;
plot([min(tx2) max(tx2)],[0 0],'k--');
m = mean(frZsc(uIx1,:),1);
se = std(frZsc(uIx1,:),0,1)./sqrt(length(uIx1)-1);
m = conv(m,gausswin(3),'same')./sum(gausswin(3));
se = conv(se,gausswin(3),'same')./sum(gausswin(3));
jbfill(tx2,m-se,m+se,[21 171 0]./255,[21 171 0]./255,1,.5);
hold on;
plot(tx2,m,'-','Color',[1 1 1]);
m = mean(frZsc(uIx2,:),1);
se = std(frZsc(uIx2,:),0,1)./sqrt(length(uIx1)-1);
m = conv(m,gausswin(3),'same')./sum(gausswin(3));
se = conv(se,gausswin(3),'same')./sum(gausswin(3));
jbfill(tx2,m-se,m+se,[255 159 0]./255,[255 159 0]./255,1,.5);
hold on;
plot(tx2,m,'-','Color',[.9 0 0]);
xlim([-500 4500]);
set(gca,'XTick',[0 2000 4000]);
xlabel('Time (ms)');
ylabel('Normalized firing rate (\sigma)');
set(gca,'LineWidth',3);
set(gca,'FontSize',14);
set(gcf,'Color','w');

%%
uIx1 = [];
uIx2 = [];
for it = 1:size(iFRzsc,1)
    
    [~,ix] = max(iFRzsc(it,:));
    if tx(ix) > 2e3 
        uIx2 = [uIx2 it];
    else
        uIx1 = [uIx1 it];
    end;
    
end;
kw = 1e3;

% uIx1 = 1:size(iFRzsc,1);
% uIx2 = [];
figure;
hold on;
plot([min(tx2) max(tx2)],[0 0],'k--');

M = mean(iFRzsc(uIx1,:),1);
SE = std(iFRzsc(uIx1,:),0,1)./sqrt(length(uIx1)-1);
M = conv(M,gausswin(kw),'same')./sum(gausswin(kw));
SE = conv(SE,gausswin(kw),'same')./sum(gausswin(kw));

jbfill(tx,M-SE,M+SE,[21 171 0]./255,[21 171 0]./255,1,.5);
hold on;
plot(tx,mean(M,1),'g');

M = mean(iFRzsc(uIx2,:),1);
SE = std(iFRzsc(uIx2,:),0,1)./sqrt(length(uIx1)-1);
M = conv(M,gausswin(kw),'same')./sum(gausswin(kw));
SE = conv(SE,gausswin(kw),'same')./sum(gausswin(kw));

jbfill(tx,M-SE,M+SE,[255 159 0]./255,[255 159 0]./255,1,.5);
hold on;
plot(tx,M,'r');
xlim([-500 4500]);
set(gca,'XTick',[0 2000 4000]);
set(gca,'YTick',[min(get(gca,'YTick')) max(get(gca,'YTick'))]);
xlabel('Time (ms)');
ylabel('Normalized firing rate (\sigma)');
set(gca,'LineWidth',3);
set(gca,'FontSize',14);
set(gcf,'Color','w');

%%
uIx1LFP1= [];uIx1LFP2= [];



cnt = 0;
for it = 1:length(uIx1)
    
    if any(strcmp(unitPoolLabel(uIx1(it)),lfpPoolLabel))
        cnt = cnt+1;
        uIx1LFP1(cnt) = uIx1(it);
        uIx1LFP2(cnt) = find(strcmp(unitPoolLabel(uIx1(it)),lfpPoolLabel));
    end;
    
end;

uIx2LFP1= [];uIx2LFP2= [];
cnt = 0;
for it = 1:length(uIx2)
    
    if any(strcmp(unitPoolLabel(uIx2(it)),lfpPoolLabel))
        cnt = cnt+1;
        uIx2LFP1(cnt) = uIx2(it);
        uIx2LFP2(cnt) = find(strcmp(unitPoolLabel(uIx2(it)),lfpPoolLabel));
    end;
end;

%%
T = 1;
W = 1/T;
TW = T*W;
k = 2*TW-1;

params1                  = [];
params1.Fs               = 1e3;
params1.pad              = 5;
params1.fpass            = [0.5 30];
params1.tapers           = [TW k];
params1.trialave         = 1;

T = 1;
W = 4;
TW = T*W;
k = 2*TW-1;

params2                  = [];
params2.Fs               = 1e3;
params2.pad              = 5;
params2.fpass            = [30 256];
params2.tapers           = [TW k];
params2.trialave         = 1;

S1b = zeros(length(LFPpoolH),967);
S1p = zeros(length(LFPpoolH),967);

S2b = zeros(length(LFPpoolH),967);
S2p = zeros(length(LFPpoolH),967);
tIx1 = find(-1000:5000 >=-1000 & -1000:5000 <=0);
tIx2 = find(-1000:5000 >=500 & -1000:5000 <=1500);
tIx3 = find(-1000:5000 >=2500 & -1000:5000 <=3500);
for jt = 1:length(LFPpoolH)
    fprintf([num2str(jt),'/',num2str(length(LFPpoolH))]);    
    x = LFPpoolH{jt};
    x = x-ones(size(x,1),1)*mean(x,1);
     
    [s1,~] = mtspectrumc((x(tIx1,:)')', params1);
    [s2,fx3] = mtspectrumc((x(tIx2,:)')', params1);
    S1p(jt,:) = s2;
    S1b(jt,:) = s1;    
   
    [s1,~] = mtspectrumc((x(tIx1,:)')', params1);
    [s2,fx4] = mtspectrumc((x(tIx3,:)')', params1);
    
    S2b(jt,:) = s1;
    S2p(jt,:) = s2;
    fprintf('\n');
end;

%%
T = 1;
W = 1/T;
TW = T*W;
k = 2*TW-1;

params1                  = [];
params1.Fs               = 1e3;
params1.pad              = 5;
params1.fpass            = [0.5 30];
params1.tapers           = [TW k];
params1.trialave         = 0;

T = 1;
W = 4;
TW = T*W;
k = 2*TW-1;

params2                  = [];
params2.Fs               = 1e3;
params2.pad              = 5;
params2.fpass            = [30 256];
params2.tapers           = [TW k];
params2.trialave         = 0;

tIx1 = find(-1000:5000 >=-1000 & -1000:5000 <=0);
tIx2 = find(-1000:5000 >=0 & -1000:5000 <=1000);
x = avgLFPh;
x = x-ones(size(x,1),1)*mean(x,1);

[s1,~] = mtspectrumc((x(tIx1,:)')', params1);
[s2,fx3] = mtspectrumc((x(tIx2,:)')', params1);
S3p = s2;
S3b = s1;

tIx1 = find(-1000:5000 >=-1000 & -1000:5000 <=0);
tIx2 = find(-1000:5000 >=2000 & -1000:5000 <=3000);
[s1,~] = mtspectrumc((x(tIx1,:)')', params1);
[s2,fx4] = mtspectrumc((x(tIx2,:)')', params1);

S4b = s1;
S4p = s2;

%%
Fs = 1e3;

%lower frequencies
movingwin1 = [1 .1];%.2
T = movingwin1(1);% length of time window in s
W = 1/T;% smoothing (+/- 1Hz)
TW = T*W;% time-bandwidth product
k = floor(2*(T*W)-1);% number of tapers

paramsTF1              = [];
paramsTF1.Fs           = Fs;
paramsTF1.pad          = 5;
paramsTF1.fpass        = [0.5 30];
paramsTF1.tapers       = [TW k];
paramsTF1.trialave     = 0;

%higher frequencies
movingwin2 = [0.25 0.0625];
T = movingwin2(1);% length of time window in s
W = 8;% smoothing (+/- 10Hz)cd
TW = T*W;% time-bandwidth product
k = floor(2*(T*W)-1);% number of tapers

paramsTF2              = [];
paramsTF2.Fs           = Fs;
paramsTF2.pad          = 5;
paramsTF2.fpass        = [30 260];
paramsTF2.tapers       = [TW k];
paramsTF2.trialave     = 0;

%%
[stPOW1] = cell(1,length(LFPpoolH));
[stPOW2] = cell(1,length(LFPpoolH));
for it = 1:length(LFPpoolH)
    fprintf([num2str(it),'/',num2str(length(LFPpoolH))]);
    x = (LFPpoolH{it}')';
    x = x-ones(size(x,1),1)*mean(x,1);
    
    [S,t1,f1] = mtspecgramc( x, movingwin1, paramsTF1);
    t1 = t1-1;
    S = S(t1>=-.4 & t1 <=5,:,:);
    t1 = t1(t1>=-.4 & t1 <=5);
    stPOW1{it} = S;
    
    [S,t2,f2] = mtspecgramc( x, movingwin2, paramsTF2);
    t2 = t2-1;
    S = S(t2>=-.4 & t2 <=5,:,:);
    t2 = t2(t2>=-.4 & t2 <=5);
    stPOW2{it} = S;
    fprintf('\n');
end;

[S,t1,f1] = mtspecgramc( (avgLFPh-ones(size(avgLFPh,1),1)*mean(avgLFPh,1) ), movingwin1, paramsTF1);
t1 = t1-1;
S = S(t1>=-.4 & t1 <=5,:,:);
t1 = t1(t1>=-.4 & t1 <=5);
erpPOW1 = S;

[S,t2,f2] = mtspecgramc( (avgLFPh-ones(size(avgLFPh,1),1)*mean(avgLFPh,1) ), movingwin2, paramsTF2);
t2 = t2-1;
S = S(t2>=-.4 & t2<=5,:,:);
t2 = t2(t2>=-.4&t2<=5);
erpPOW2 = S;

% %%
% stPOW3 = cell(1,length(LFPpoolM));
% % stPOW4 = cell(1,length(LFPpoolM));
% 
% for it = 1:length(LFPpoolM)
%     
%     x = LFPpoolM{it};
%     x = x-ones(size(x,1),1)*mean(x,1);
%     
%     [S,t1,f1] = mtspecgramc( x, movingwin1, paramsTF1);
%     t1 = t1-1;
%     S = S(t1>=-.4 & t1 <=4.5,:,:);
%     t1 = t1(t1>=-.4 & t1 <=4.5);
%     stPOW3{it} = S;
%     
% %     [S,t2,f2] = mtspecgramc( x, movingwin2, paramsTF2);
% %     t2 = t2-1;
% %     S = S(t2>=-.4 & t2 <=4.5,:,:);
% %     t2 = t2(t2>=-.4 & t2 <=4.5);
% %     stPOW4{it} = S;
%     
% end;
% 
% 
% % [S,t1,f1] = mtspecgramc( (avgLFPm-ones(size(avgLFPm,1),1)*mean(avgLFPm,1) ), movingwin1, paramsTF1);
% % t1 = t1-1;
% % S = S(t1>=-.4 & t1 <=5,:,:);
% % t1 = t1(t1>=-.4 & t1 <=5);
% % erpPOW3 = S;
% % 
% % 
% % [S,t2,f2] = mtspecgramc( gradient(avgLFPm')', movingwin2, paramsTF2);
% % t2 = t2-1;
% % S = S(t2>=-.75 & t2<=4.5,:,:);
% % t2 = t2(t2>=-.75&t2<=4.5);
% % erpPOW4 = S;
% 
% %%
% for it = 1:length(stPOW1)
%     figure;
%     imagesc(t1,f1,20*log10(squeeze(mean(stPOW1{(it)},3)))');
%     axis xy;colorbar;
%     colormap jet;
% end;
% 
% %%
% for it = 1:size(X2,1)
%     figure;
%     imagesc(t1,f1,squeeze(log10(X2(it,:,:))-log10(M2(it,:,:)))');
%     axis xy;colorbar;
%     colormap jet;
% end;

%%
Fs = 1e3;
[nsmp,~] = size(LFPpoolH{1});

cfg                     = [];
cfg.foi                 = 1.5:0.25:30;%;
cfg.toi                 = [0:nsmp*3-1]./Fs;
cfg.pad                 = 'nextpow2'; 
cfg.method              = 'wavelet';
cfg.width               = 4;
cfg.output              = 'fourier';

%sigSelItc = zeros(1,length( LFPpoolH ));
itc = zeros(length(cfg.toi)/3,length(cfg.foi),length( LFPpoolH ));
parfor it = 1:length( LFPpoolH )
    
    [nsmp,ntrl] = size(LFPpoolH{it});
    
    dum = [];
    dum.label = {'dumChan1'};
    dum.trial = cell(1,ntrl);
    dum.time = cell(1,ntrl);
    for jt = 1:ntrl
        x = LFPpoolH{it}(:,jt)';
        dum.trial{jt} = [fliplr(x) x fliplr(x)];
        dum.time{jt} = cfg.toi;
    end;
    
    [phi] = ft_freqanalysis( cfg , dum);        
    phi.fourierspctrm = phi.fourierspctrm(:,:,:,nsmp+1:nsmp*2);
    phi.time = phi.time(nsmp+1:nsmp*2);
    phi.time = phi.time-7;
    
%     tmp = angle(squeeze(phi.fourierspctrm));
%     concat = zeros(size(tmp,1),size(tmp,2)*size(tmp,3));
%     ix = 1:size(tmp,3);
%     for kt = 1:size(tmp,2)
%         concat(:,ix) = squeeze(tmp(:,kt,:));       
%         ix = ix+size(tmp,3);
%     end;
%     p = zeros(1,size(concat,2));
%     z = zeros(1,size(concat,2));
%     for kt = 1:size(concat,2)
%         [p(kt),z(kt)] = circ_rtest(concat(:,kt));
%     end;
%     p = reshape(p,[length(phi.time) length(phi.freq)]);
%     z = reshape(z,[length(phi.time) length(phi.freq)]);
%     
%     [hcItc]=fdr_bh(p,0.001);
%     if any(hcItc==1)
%         sigSelItc(it) = 1;
%     end;
    
    tmp = phi.fourierspctrm;
    tmp = tmp./abs(tmp);
    tmp = sum(tmp,1);
    tmp = abs(tmp)/ntrl;
    tmp = squeeze(tmp);
    itc(:,:,it) = tmp';
    
end;
itcTime = cfg.toi(nsmp+1:nsmp*2)-7;
itcFreq = cfg.foi;

%%
X = zeros(length(stPOW1),length(t1),length(f1));
M = zeros(length(stPOW1),length(t1),length(f1));
for it = 1:length(stPOW1)    
    x = (squeeze(mean(stPOW1{(it)},3)));
    m = ones(size(x,1),1)*mean(x(t1<=0,:),1);
    X(it,:,:) = x;
    M(it,:,:) = m;    
end;

% X2 = zeros(length(stPOW3),length(t1),length(f1));
% M2 = zeros(length(stPOW3),length(t1),length(f1));
% for it = 1:length(stPOW3)    
%     x = (squeeze(mean(stPOW3{(it)},3)));
%     m = ones(size(x,1),1)*mean(x(t1<=0,:),1);
%     X2(it,:,:) = x;
%     M2(it,:,:) = m;    
% end;

empT = zeros(length(t1),length(f1));
pvalT = zeros(length(t1),length(f1));

for it = 1:length(t1)
    for jt = 1:length(f1)
        
        d = X(:,it,jt)-M(:,it,jt);
        [~,p,~,stat] = ttest(d);
        empT(it,jt) = stat.tstat;
        pvalT(it,jt) = p;
    end;
end;
[hcT]=fdr_bh(pvalT,0.01);

%%
dumIx1 = [uIx1LFP1 uIx2LFP1];
dumIx2 = [uIx1LFP2 uIx2LFP2];

a = zeros(1,size(iFRzsc,1));
b = zeros(1,size(iFRzsc,1));
for it = 1:size(iFRzsc,1)
    
    a(it) = max(iFRzsc(it,tx>=0 & tx<2e3));
    b(it) = max(iFRzsc(it,tx>=2e3 & tx<=4e3));
       
end;

fr1 = zeros(1,length(dumIx1));
fr2 = zeros(1,length(dumIx1));
p1 = zeros(1,length(dumIx1));
p2 = zeros(1,length(dumIx1));
p3 = zeros(1,length(dumIx1));
p4 = zeros(1,length(dumIx1));
for it = 1:length(dumIx1)
    
    fr1(it) = max(iFRzsc(dumIx1(it),tx>=0 & tx<2e3));
    fr2(it) = max(iFRzsc(dumIx1(it),tx>=2e3 & tx<=4e3));
    
    p1(it) = squeeze(mean(mean(X(dumIx2(it),t1>=0 & t1<1,f1>=1 & f1<=6),3),2));
    p2(it) = squeeze(mean(mean(X(dumIx2(it),t1>=2 & t1<=3,f1>=1 & f1<=6),3),2));
    
    p3(it) = squeeze(mean(mean(erpPOW1(t1>=0 & t1<2,f1>=1 & f1<=6,dumIx2(it)),2),1));
    p4(it) = squeeze(mean(mean(erpPOW1(t1>=2 & t1<=4,f1>=3 & f1<=6,dumIx2(it)),2),1));
    
end;

%%
figure;
ax = gca;
hold on;
x = normalize(a);
y = normalize(b);
plot(x(uIx1),y(uIx1),'o','Color',[21 171 0]./255);
plot(x(uIx2),y(uIx2),'x','Color',[255 159 0]./255);
plot([0 1],[0 1],'k');
rxy1 = corr(x(uIx1)',y(uIx1)','Type','Spearman');
rxy2 = corr(x(uIx2)',y(uIx2)','Type','Spearman');
%title(['rxy1: ',num2str(round(rxy1*100)/100),' rxy2: ',num2str(round(rxy2*100)/100)]);

figure;
ax = [ax gca];
hold on;
for it = 1:length(uIx1)
    plot([1 2],[x(uIx1(it)) y(uIx1(it))],'Color',[.75 .75 .75]);    
end;
for it = 1:length(uIx2)
    plot([3 4],[x(uIx2(it)) y(uIx2(it))],'Color',[.75 .75 .75]);
end;
plot(ones(1,length(uIx1)),x(uIx1),'o','Color',[21 171 0]./255);
plot(2*ones(1,length(uIx1)),y(uIx1),'o','Color',[21 171 0]./255);
plot(3*ones(1,length(uIx2)),x(uIx2),'o','Color',[255 159 0]./255);
plot(4*ones(1,length(uIx2)),y(uIx2),'o','Color',[255 159 0]./255);
xlim([0 5]);

figure;
ax = [ax gca];
hold on;
x = normalize(fr1);
y = normalize(p1);
plot(x(ismember(dumIx1,uIx1LFP1)),y(ismember(dumIx1,uIx1LFP1)),'o','Color',[21 171 0]./255);
plot(x(ismember(dumIx1,uIx2LFP1)),y(ismember(dumIx1,uIx2LFP1)),'x','Color',[255 159 0]./255);
plot([0 1],[0 1],'k','LineSmoothing','on');
rxy1 = corr(x(ismember(dumIx1,uIx1LFP1))',y(ismember(dumIx1,uIx1LFP1))','Type','Spearman');
rxy2 = corr(x(ismember(dumIx1,uIx2LFP1))',y(ismember(dumIx1,uIx2LFP1))','Type','Spearman');
%title(['rxy1: ',num2str(round(rxy1*100)/100),' rxy2: ',num2str(round(rxy2*100)/100)]);

figure;
ax = [ax gca];
hold on;
x = normalize(fr2);
y = normalize(p2);
plot(x(ismember(dumIx1,uIx1LFP1)),y(ismember(dumIx1,uIx1LFP1)),'o','Color',[21 171 0]./255);
plot(x(ismember(dumIx1,uIx2LFP1)),y(ismember(dumIx1,uIx2LFP1)),'x','Color',[255 159 0]./255);
plot([0 1],[0 1],'k','LineSmoothing','on');
rxy1 = corr(x(ismember(dumIx1,uIx1LFP1))',y(ismember(dumIx1,uIx1LFP1))','Type','Spearman');
rxy2 = corr(x(ismember(dumIx1,uIx2LFP1))',y(ismember(dumIx1,uIx2LFP1))','Type','Spearman');
%title(['rxy1: ',num2str(round(rxy1*100)/100),' rxy2: ',num2str(round(rxy2*100)/100)]);

figure;
ax = [ax gca];
hold on;
x = normalize(fr1);
y = normalize(p3);
plot(x(ismember(dumIx1,uIx1LFP1)),y(ismember(dumIx1,uIx1LFP1)),'o','Color',[21 171 0]./255);
plot(x(ismember(dumIx1,uIx2LFP1)),y(ismember(dumIx1,uIx2LFP1)),'x','Color',[255 159 0]./255);
plot([0 1],[0 1],'k','LineSmoothing','on');
rxy1 = corr(x(ismember(dumIx1,uIx1LFP1))',y(ismember(dumIx1,uIx1LFP1))','Type','Spearman');
rxy2 = corr(x(ismember(dumIx1,uIx2LFP1))',y(ismember(dumIx1,uIx2LFP1))','Type','Spearman');
%title(['rxy1: ',num2str(round(rxy1*100)/100),' rxy2: ',num2str(round(rxy2*100)/100)]);

figure;
ax = [ax gca];
hold on;
x = normalize(fr2);
y = normalize(p4);
plot(x(ismember(dumIx1,uIx1LFP1)),y(ismember(dumIx1,uIx1LFP1)),'o','Color',[21 171 0]./255);
plot(x(ismember(dumIx1,uIx2LFP1)),y(ismember(dumIx1,uIx2LFP1)),'x','Color',[255 159 0]./255);
plot([0 1],[0 1],'k','LineSmoothing','on');
rxy1 = corr(x(ismember(dumIx1,uIx1LFP1))',y(ismember(dumIx1,uIx1LFP1))','Type','Spearman');
rxy2 = corr(x(ismember(dumIx1,uIx2LFP1))',y(ismember(dumIx1,uIx2LFP1))','Type','Spearman');
%title(['rxy1: ',num2str(round(rxy1*100)/100),' rxy2: ',num2str(round(rxy2*100)/100)]);

set(ax,'LineWidth',2);
set(ax,'Fontsize',14);
for it = 1:length(ax)
    set(ax(it),'YTick',[min(get(ax(it),'YTick')) max(get(ax(it),'YTick'))]);
    set(ax(it),'XTick',[min(get(ax(it),'XTick')) max(get(ax(it),'XTick'))]);
    box(ax(it),'on');
    xlabel(ax(it),'Firing Rate (a.u.)');
    if ismember(it,1:2)
        ylabel(ax(it),'Firing Rate (a.u.)');
    else
        ylabel(ax(it),'LFP-power (a.u.)');
    end;
end;
set(ax(2),'XTick',1:4);
set(ax(2),'XTickLabel',[]);

savepath = '/home/rouxf/tmp/figuresLNM/';
for it = 1:get(gcf,'Number')
    set(it,'Color','w');
    saveas(it,[savepath,'PowerVsFiringRate',num2str(it),'.fig']);
end;


%%
figure;
a = gca;
hold on;
imagesc(t1.*1e3,f1,(squeeze(mean((10*log10(X)-10*log10(M)),1)).*(hcT))');
colormap jet;axis xy;
plot([0 0],[min(f1) max(f1)],'k');
plot([2 2].*1e3,[min(f1) max(f1)],'k');
axis tight;
cb = colorbar;
set(get(cb,'YLabel'),'String',{'LFP-power','(dB)'});
set(cb,'YTick',round([min(get(cb,'YTick')) max(get(cb,'YTick'))].*10)/10);
caxis([round([min(get(cb,'YTick')) max(get(cb,'YTick'))].*10)/10]);

figure;
a = [a gca];
ERPm = zeros(size(erpPOW1));
for it = 1:size(erpPOW1,3)
    ERPm(:,:,it) = ones(size(erpPOW1(:,:,it),1),1)*squeeze(mean(erpPOW1(t1<=0,:,it),1));
end;
hold on;
imagesc(t1.*1e3,f1,squeeze(mean((10*log10(erpPOW1)-10*log10(ERPm)),3))');
plot([0 0],[min(f1) max(f1)],'k');
plot([2 2].*1e3,[min(f1) max(f1)],'k');
axis xy;axis tight;
cb = colorbar;
set(get(cb,'YLabel'),'String',{'LFP-power','(dB)'});
set(cb,'YTick',round([min(get(cb,'YTick')) max(get(cb,'YTick'))].*10)/10);
caxis([round([min(get(cb,'YTick')) max(get(cb,'YTick'))].*10)/10]);

figure;
a = [a gca];
hold on;
imagesc(t1.*1e3,f1,(empT.*(hcT))');%pval<0.05/length(empT(:)))
plot([0 0],[min(f1) max(f1)],'k');
plot([2 2].*1e3,[min(f1) max(f1)],'k');
axis xy;caxis([-max(max(abs(empT))) max(max(abs(empT)))]);
colormap jet;axis tight;
cb = colorbar;
set(get(cb,'YLabel'),'String',{'T-values','(p<0.01, corrected)'});
set(cb,'YTick',round([min(get(cb,'YTick')) max(get(cb,'YTick'))].*10)/10);
caxis([round([min(get(cb,'YTick')) max(get(cb,'YTick'))].*10)/10]);

figure;
a = [a gca];
hold on;
imagesc(itcTime.*1e3,itcFreq,squeeze(mean(itc,3))');
shading interp;
plot([0 0],[min(itcFreq) max(itcFreq)],'k');
plot([2e3 2e3],[min(itcFreq) max(itcFreq)],'k');
axis tight;xlim([min(t1) max(t1)].*1e3);
cb = colorbar;
set(get(cb,'YLabel'),'String',{'LFP','Inter-trial coherence (a.u.)'});
set(cb,'YTick',round([min(get(cb,'YTick')) max(get(cb,'YTick'))].*10)/10);
caxis([round([min(get(cb,'YTick')) max(get(cb,'YTick'))].*10)/10]);

figure;
a = [a gca];
x = 10*log10(erpPOW1)-10*log10(ERPm);
m = mean(mean(x(:,f1>=3 & f1<=5,:),2),3)';
se = (std(mean(x(:,f1>=3 & f1<=5,:),2),0,3)./sqrt(size(x,3)-1))';
plot([min(t1) max(t1)].*1e3,[0 0],'k--');
jbfill(t1.*1e3,m-se,m+se,'r','r',1,.5);
hold on;
plot(t1.*1e3,m,'r');
m = squeeze(mean(mean(10*log10(X(:,:,f1>=3 & f1<=5)),3)-mean(10*log10(M(:,:,f1>=3 & f1<=5)),3),1));
se = squeeze(std(mean(10*log10(X(:,:,f1>=3 & f1<=5)),3)-mean(10*log10(M(:,:,f1>=3 & f1<=5)),3),0,1))./sqrt(size(X,1)-1);
jbfill(t1.*1e3,m-se,m+se,'c','c',1,.5);
hold on;
plot(t1.*1e3,m,'b');
axis tight;

figure;
a = [a gca];
x = 10*log10(S1p);
m = mean(x,1);        
se = std(x,0,1)./sqrt(size(x,1)-1);
jbfill((fx3),m-se,m+se,'c','c',1,.5);
hold on;
plot((fx3),m,'Color','b');
x = 10*log10(S2p);
m = mean(x,1);        
se = std(x,0,1)./sqrt(size(x,1)-1);
jbfill((fx3),m-se,m+se,'r','r',1,.5);
hold on;
plot((fx3),m,'Color','r');
x = 10*log10(S1b);
m = mean(x,1);        
se = std(x,0,1)./sqrt(size(x,1)-1);
jbfill((fx3),m-se,m+se,'k','k',1,.5);
hold on;
plot((fx3),m,'Color','k');
axis tight;

figure;
a = [a gca];
x = (S1p);
m = mean(x,1);
y = log10(m)';
b = regress(y,[ones(size(fx3')) log10(fx3)']);
yp = [ones(size(fx3')) log10(fx3)']*b;
m = 10.^(y-yp);
se = std(x,0,1)./sqrt(size(x,1)-1);
y = log10(se)';
b = regress(y,[ones(size(fx3')) log10(fx3)']);
yp = [ones(size(fx3')) log10(fx3)']*b;
se = 10.^(y-yp);
%jbfill(fx3,m'-se',m'+se','c','c',1,.5);hold on;
plot((fx3),m,'Color','b');
x = (S2p);
m = mean(x,1);
y = log10(m)';
b = regress(y,[ones(size(fx3')) log10(fx3)']);
yp = [ones(size(fx3')) log10(fx3)']*b;
m = 10.^(y-yp);
se = std(x,0,1)./sqrt(size(x,1)-1);
y = log10(se)';
b = regress(y,[ones(size(fx3')) log10(fx3)']);
yp = [ones(size(fx3')) log10(fx3)']*b;
se = 10.^(y-yp);
hold on;
%jbfill(fx3,m'-se',m'+se','r','r',1,.5);hold on;
plot((fx3),m,'Color','r');
x = (S1b);
m = mean(x,1);
y = log10(m)';
b = regress(y,[ones(size(fx3')) log10(fx3)']);
yp = [ones(size(fx3')) log10(fx3)']*b;
m = 10.^(y-yp);
plot((fx3),m,'Color','k');
plot([min(fx3) max(fx3)],[0 0],'k--');
axis tight;

set(a,'LineWidth',3);
set(a,'Fontsize',14);
set(gcf,'Color','w');
set(a([1 2 3 4 5]),'XTick',[0 2 4].*1e3);
for it = [1 2 3 4 5]
    xlabel(a(it),'Time (ms)');
    ylabel(a(it),'Frequency (Hz)');
    if ismember(it,1:4)
        set(a(it),'YTick',[5 15 25]);
    end;
end;
ylabel(a(5),{'Power (dB)';'LFP'});
xlabel(a(6),'Frequency (Hz)');
ylabel(a(6),{'Power (log)';'LFP'});
xlabel(a(7),'Frequency (Hz)');
ylabel(a(7),{'Power (1/f corrected)';'LFP'});

for it = [1 2 3 4]
    box(a(it),'on');
end;
for it = [5 6 7]
    box(a(it),'off');
    set(a(it),'YTick',[min(get(a(it),'YTick')) max(get(a(it),'YTick'))]);
    ylim(a(it),[min(get(a(it),'YTick')) max(get(a(it),'YTick'))]);
end;

for curFig = 1:get(gcf,'Number');
    set(curFig,'Color','w');
    colormap(curFig,'jet');
end;

savepath = '/home/rouxf/tmp/figuresLNM/';
for it = 1:get(gcf,'Number')
    set(it,'Color','w');
    saveas(it,[savepath,'LFPpowerAnalysis',num2str(it),'.fig']);
end;

% %%
% X2 = zeros( length(stPOW2), length(t2),length(f2) );
% M2 = zeros( length(stPOW2), length(t2),length(f2) );
% 
% for it = 1:length(stPOW2)
%     x = squeeze(mean( stPOW2{(it)} ,3));    
%     m = ones(size(x,1),1)*mean(x(t2 <=0,:),1);
%     M2(it,:,:) = m;
%     X2(it,:,:) = x;
% end;
% 
% figure;
% imagesc(t2,f2,squeeze(mean(10*log10(X2)-10*log10(M2),1))');
% axis xy;
% colormap jet;
% 
% pvalT2 = zeros(length(t2),length(f2));
% empT2 = zeros(length(t2),length(f2));
% for it = 1:length(t2)
%     for jt = 1:length(f2)
%         d = squeeze(X2(:,it,jt))-squeeze(M2(:,it,jt));
%         [~,p,~,stat] = ttest(d);
%         pvalT2(it,jt) = p;
%         empT2(it,jt) = stat.tstat;
%         
%     end;
% end;
% 
% [hcT2, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvalT2,0.05);
% 
% figure;
% subplot(121);
% imagesc(t2,f2,(empT2'));
% axis xy;
% subplot(122);
% imagesc(t2,f2,(empT2.*(hcT2))');
% axis xy;

%%
Fs = 1e3;
[nsmp,~] = size(LFPpoolH{1});

cfg                     = [];
cfg.foi                 = 2.^((0:40)./8);
cfg.toi                 = [0:nsmp*3-1]./Fs;
cfg.pad                 = 'nextpow2'; 
cfg.method              = 'wavelet';
cfg.width               = 4;
cfg.output              = 'fourier';

phi = cell(1,length(LFPpoolH));
parfor curMW = 1:length( LFPpoolH )
    
    [nsmp,ntrl] = size(LFPpoolH{curMW});
    
    dum = [];
    dum.label = {'dumChan1'};
    dum.trial = cell(1,ntrl);
    dum.time = cell(1,ntrl);
    for jt = 1:ntrl
        x = LFPpoolH{curMW}(:,jt)';
        dum.trial{jt} =[fliplr(x) x fliplr(x)];
        dum.time{jt} = cfg.toi;
    end;
        
    [phi{curMW}] = ft_freqanalysis( cfg , dum);        
    phi{curMW}.fourierspctrm = phi{curMW}.fourierspctrm(:,:,:,nsmp+1:nsmp*2);
    phi{curMW}.time = phi{curMW}.time(nsmp+1:nsmp*2);
    phi{curMW}.time = -1:1e-3:5;           
end;

%%
for it = 1:2:size(itc,3)
figure;
subplot(421);
a = gca;
hold on;
pcolor(itcTime.*1e3,itcFreq,squeeze(itc(:,:,it))');
ix = find(strcmp(lfpPoolLabel(it),unitPoolLabel));
plot(tx2,normalize(frPool(ix(1),:)).*10+2,'w','LineWidth',2);
plot([0 0],[min(itcFreq) max(itcFreq)],'k');
plot([2 2].*1e3,[min(itcFreq) max(itcFreq)],'k');
axis tight;xlim([-400 4.5e3]);
shading interp;

subplot(422);
a = [a gca];
hold on;
pcolor(itcTime.*1e3,itcFreq,squeeze(itc(:,:,it+1))');
ix = find(strcmp(lfpPoolLabel(it+1),unitPoolLabel));
plot(tx2,normalize(frPool(ix(1),:)).*10+2,'w','LineWidth',2);
plot([0 0],[min(itcFreq) max(itcFreq)],'k');
plot([2 2].*1e3,[min(itcFreq) max(itcFreq)],'k');
axis tight;xlim([-400 4.5e3]);
shading interp;
xlim([-400 4.5e3]);
colormap jet;

subplot(423);
a = [a gca];
x = LFPpoolH{it};
z = mean(sqrt(x.^2),1);
z = (z-mean(z))./std(z);
x(:,abs(z)>2) = [];
bpf = [2.5 9];
b = fir1(3*fix(1e3/bpf(1)),bpf./(1e3/2),'bandpass');
f = filtfilt(b,1,mean(x,2));
hold on;
plot(-1e3:5e3,normalize(mean(x,2)),'k');
plot(-1e3:5e3,normalize(f),'g','LineWidth',1);
xlim([-400 4.5e3]);
axis off;

subplot(424);
a = [a gca];
x = LFPpoolH{it+1};
z = mean(sqrt(x.^2),1);
z = (z-mean(z))./std(z);
x(:,abs(z)>2) = [];
bpf = [2.5 9];
b = fir1(3*fix(1e3/bpf(1)),bpf./(1e3/2),'bandpass');
f = filtfilt(b,1,mean(x,2));
hold on;
plot(-1e3:5e3,normalize(mean(x,2)),'k');
plot(-1e3:5e3,normalize(f),'g','LineWidth',1);
xlim([-400 4.5e3]);
axis off;

subplot(425);
a = [a gca];
hold on;
ix = find(strcmp(lfpPoolLabel(it),unitPoolLabel));
spk = spkPool{ix(1)};
for jt = 1:size(spk,1)
    ix = spk(jt,:)==1;
    ts = tx(ix);
    y = jt*ones(1,length(ts));
    ts = [ts;ts];
    y = [y-.5;y+.5];
    h = line(ts,y,'Color','k');
    set(h,'LineWidth',.1);
end;
plot([0 0],[1 size(spk,1)],'r');
plot([2 2].*1e3,[1 size(spk,1)],'r');
ylim([1 size(spk,1)]);
xlim([-400 4.5e3]);

subplot(426);
a = [a gca];
hold on;
ix = find(strcmp(lfpPoolLabel(it+1),unitPoolLabel));
spk = spkPool{ix(1)};
for jt = 1:size(spk,1)
    ix = spk(jt,:)==1;
    ts = tx(ix);
    y = jt*ones(1,length(ts));
    ts = [ts;ts];
    y = [y-.5;y+.5];
    h = line(ts,y,'Color','k');
    set(h,'LineWidth',.1);
end;
plot([0 0],[1 size(spk,1)],'r');
plot([2 2].*1e3,[1 size(spk,1)],'r');
ylim([1 size(spk,1)]);
xlim([-400 4.5e3]);

subplot(427); 
a = [a gca];
x = LFPpoolH{it};
z = mean(sqrt(x.^2),1);
z = (z-mean(z))./std(z);
x(:,abs(z)>2) = [];
ix = find(strcmp(lfpPoolLabel(it),unitPoolLabel));
spk = spkPool{ix(1)};
if (delIxPool{it}~=0)
    spk(delIxPool{it},:) = [];
end;
spk = spk(ixHPool{it},:)';
spk(:,abs(z)>2) = [];
lag = 500;
cnt1 = 0;
sta1 = [];
for jt = 1:size(spk,2)
    tsIx = find(spk(:,jt));
    for kt = 1:length(tsIx)
        if (tx(tsIx(kt))-lag>0) && (tx(tsIx(kt))+lag<=2e3)
            cnt1 = cnt1+1;
            sta1(cnt1,:) = x(tsIx(kt)-lag:tsIx(kt)+lag,jt);
        end;
    end;
end;
cnt2 = 0;
sta2 = [];
for jt = 1:size(spk,2)
    tsIx = find(spk(:,jt));
    for kt = 1:length(tsIx)
        if (tx(tsIx(kt))-lag>2e3) && (tx(tsIx(kt))+lag<=4e3)
            cnt2 = cnt2+1;
            sta2(cnt2,:) = x(tsIx(kt)-lag:tsIx(kt)+lag,jt);
        end;
    end;
end;
hold on;
plot(-lag:lag,mean(sta1,1),'b');
plot(-lag:lag,mean(sta2,1),'r');
axis tight;

subplot(428); 
a = [a gca];
x = LFPpoolH{it+1};
z = mean(sqrt(x.^2),1);
z = (z-mean(z))./std(z);
x(:,abs(z)>2) = [];
ix = find(strcmp(lfpPoolLabel(it+1),unitPoolLabel));
spk = spkPool{ix(1)};
if (delIxPool{it+1}~=0)
    spk(delIxPool{it+1},:) = [];
end;
spk = spk(ixHPool{it+1},:)';
spk(:,abs(z)>2) = [];
lag = 500;
cnt1 = 0;
sta1 = [];
for jt = 1:size(spk,2)
    tsIx = find(spk(:,jt));
    for kt = 1:length(tsIx)
        if (tx(tsIx(kt))>0) && (tx(tsIx(kt))<=2e3)
            cnt1 = cnt1+1;
            sta1(cnt1,:) = x(tsIx(kt)-lag:tsIx(kt)+lag,jt);
        end;
    end;
end;
cnt2 = 0;
sta2 = [];
for jt = 1:size(spk,2)
    tsIx = find(spk(:,jt));
    for kt = 1:length(tsIx)
        if (tx(tsIx(kt))-lag>2e3) && (tx(tsIx(kt))<=4e3)
            cnt2 = cnt2+1;
            sta2(cnt2,:) = x(tsIx(kt)-lag:tsIx(kt)+lag,jt);
        end;
    end;
end;
hold on;
plot(-lag:lag,mean(sta1,1),'b');
plot(-lag:lag,mean(sta2,1),'r');
axis tight;

set(a(1:6),'XTick',[0 2 4].*1e3);
for jt = 1:length(a)
    xlabel(a(jt),'Time (ms)');
    if ismember(jt,[1 2])
        ylabel(a(jt),'Frequency (Hz)');
    end;
end;
set(gcf,'Color','w');

end;
savepath = '/home/rouxf/tmp/figuresLNM/';
for it = 1:get(gcf,'Number')
    set(it,'Color','w');
    saveas(it,[savepath,'lfpITCfiringRate',num2str(it),'.fig']);
end;

%%

ZvalR1 = cell(1,length( LFPpoolH));
PvalR1 = cell(1,length( LFPpoolH));
ZvalR2 = cell(1,length( LFPpoolH));
PvalR2 = cell(1,length( LFPpoolH));

tIx1 = find(tx>0 & tx <=2e3);
tIx2 = find(tx>2e3 & tx <=4e3);

cnt = 0;
unitSelIx = [];
for curClu = 1:length(unitPoolLabel)
    
    [lfpIx] = find(strcmp(lfpPoolLabel,unitPoolLabel(curClu)));
    
    if ~isempty(lfpIx)
        unitSelIx  = [unitSelIx curClu];
        cnt = cnt+1;
        tmp = angle(squeeze(phi{lfpIx}.fourierspctrm));
        concat = zeros(size(tmp,2),length(tIx1)*size(tmp,1));
        ix = 1:length(tIx1);
        for jt = 1:size(tmp,1)
            concat(:,ix) = squeeze(tmp(jt,:,tIx1));
            ix = ix+length(tIx1);
        end;
        
        spk = spkPool{curClu};
        if (delIxPool{lfpIx}~=0)
            spk(delIxPool{lfpIx},:) = [];
        end;
        spk = spk(ixHPool{lfpIx},tIx1)';
        [nsmp2,ntrl2] = size(spk);
        ts = spk(:);
        p = zeros(size(concat,1),1);
        z = zeros(size(concat,1),1);
        for ft = 1:size(concat,1)
            [p(ft),z(ft)] = circ_rtest(concat(ft,ts==1)');
        end;
        
        ZvalR1{cnt} = z;
        PvalR1{cnt} = p;
        
        tmp = angle(squeeze(phi{lfpIx}.fourierspctrm));
        concat = zeros(size(tmp,2),length(tIx2)*size(tmp,1));
        ix = 1:length(tIx2);
        for jt = 1:size(tmp,1)
            concat(:,ix) = squeeze(tmp(jt,:,tIx2));
            ix = ix+length(tIx2);
        end;
        
        spk = spkPool{curClu};
        if (delIxPool{lfpIx}~=0)
            spk(delIxPool{lfpIx},:) = [];
        end;
        spk = spk(ixHPool{lfpIx},tIx2)';
        [nsmp2,ntrl2] = size(spk);
        ts = spk(:);
        p = zeros(size(concat,1),1);
        z = zeros(size(concat,1),1);
        for ft = 1:size(concat,1)
            [p(ft),z(ft)] = circ_rtest(concat(ft,ts==1)');
        end;
        
        ZvalR2{cnt} = z;
        PvalR2{cnt} = p;
    end;
end;

cnt = 0;
PLVsigIx1 = [];
for it = 1:length(PvalR1)
    
    [hcPLV] = fdr_bh(PvalR1{it},0.001);
    %if any(any(PvalR{it} < 0.001/length(cfg.foi)))
    if any(any(hcPLV) == 1)
        cnt = cnt+1;
        PLVsigIx1(cnt) = it;
    end;
end;

cnt = 0;
PLVsigIx2 = [];
for it = 1:length(PvalR2)
    
    [hcPLV] = fdr_bh(PvalR2{it},0.001);
    %if any(any(PvalR{it} < 0.001/length(cfg.foi)))
    if any(any(hcPLV) == 1)
        cnt = cnt+1;
        PLVsigIx2(cnt) = it;
    end;
end;

avg = zeros(length(unitSelIx),length(phi{1}.freq))';
for it = 1:length(ZvalR1)    
    avg(:,it) = ZvalR1{it};
end;

%%
figure;
hold on;
m = mean(avg(:,ismember(unitSelIx,uIx1)),2)';
se = (std(avg(:,ismember(unitSelIx,uIx1)),0,2)./sqrt(sum(ismember(unitSelIx,uIx1))-1))';
jbfill(log(phi{1}.freq),m-se,m+se,[21 171 0]./255,[21 171 0]./255,1,.5);
hold on;
plot(log(phi{1}.freq),m,'Color','g');
m = mean(avg(:,ismember(unitSelIx,uIx2)),2)';
se = (std(avg(:,ismember(unitSelIx,uIx2)),0,2)./sqrt(sum(ismember(unitSelIx,uIx2))-1))';
jbfill(log(phi{1}.freq),m-se,m+se,[255 159 0]./255,[255 159 0]./255,1,.5);
hold on;
plot(log(phi{1}.freq),m,'Color','r');
set(gca,'XTick',log([2 5:5:30]));
set(gca,'XTickLabel',{'2' '5' '10' ' ' '20' ' ' '30'});
xlabel('Frequency (Hz)');
ylabel('z-statistic');
set(gca,'Fontsize',14);
set(gca,'LineWidth',3);
set(gca,'YTick',[min(get(gca,'YTick')) max(get(gca,'YTick'))]);
xlim(log([min(phi{1}.freq) max(phi{1}.freq)]));
set(gcf,'Color','w');

avg = zeros(length(unitSelIx),length(phi{1}.freq))';
for it = 1:length(ZvalR2)    
    avg(:,it) = ZvalR2{it};
end;

figure;
hold on;
m = mean(avg(:,ismember(unitSelIx,uIx1)),2)';
se = (std(avg(:,ismember(unitSelIx,uIx1)),0,2)./sqrt(sum(ismember(unitSelIx,uIx1))-1))';
jbfill(log(phi{1}.freq),m-se,m+se,[21 171 0]./255,[21 171 0]./255,1,.5);
hold on;
plot(log(phi{1}.freq),m,'Color','g');
m = mean(avg(:,ismember(unitSelIx,uIx2)),2)';
se = (std(avg(:,ismember(unitSelIx,uIx2)),0,2)./sqrt(sum(ismember(unitSelIx,uIx2))-1))';
jbfill(log(phi{1}.freq),m-se,m+se,[255 159 0]./255,[255 159 0]./255,1,.5);
hold on;
plot(log(phi{1}.freq),m,'Color','r');
set(gca,'XTick',log([2 5:5:30]));
set(gca,'XTickLabel',{'2' '5' '10' ' ' '20' ' ' '30'});
xlabel('Frequency (Hz)');
ylabel('z-statistic');
set(gca,'Fontsize',14);
set(gca,'LineWidth',3);
set(gca,'YTick',[min(get(gca,'YTick')) max(get(gca,'YTick'))]);
xlim(log([min(phi{1}.freq) max(phi{1}.freq)]));
set(gcf,'Color','w');

savepath = '/home/rouxf/tmp/figuresLNM/';
for it = 1:get(gcf,'Number')
    set(it,'Color','w');
    saveas(it,[savepath,'phaseLockingZvalue',num2str(it),'.fig']);
end;

%%
lag = 500;

T = 2;
W = 1/T;
TW = T*W;
k = 2*TW-1;

params1                  = [];
params1.Fs               = 1e3;
params1.pad              = 5;
params1.fpass            = [0.5 30];
params1.tapers           = [TW k];
params1.trialave         = 1;

T = 2;
W = 4;
TW = T*W;
k = 2*TW-1;

params2                  = [];
params2.Fs               = 1e3;
params2.pad              = 5;
params2.fpass            = [30 256];
params2.tapers           = [TW k];
params2.trialave         = 1;

R1 = zeros(length(spkPool),1);
R2 = zeros(length(spkPool),1);
S1 = zeros(length(spkPool),1934);
S2 = zeros(length(spkPool),1934);
XC1 = zeros(length(spkPool),lag*2+1);
XC2 = zeros(length(spkPool),lag*2+1);
isi1 = zeros(length(spkPool),502);
isi2 = zeros(length(spkPool),502);
dt =-1000:4500;
parfor it = 1:length(spkPool)
    fprintf([num2str(it),'/',num2str(length(spkPool))]);
    spk = spkPool{it};
    
    x1 = spk(:,tx >=0 & tx <=2e3);
    x2 = spk(:,tx >=2e3 & tx <=4e3);
    
    xt = repmat(dt',[size(x1,1) 1]);
    xt = xt(:);
    
    ix = x1';
    sts = xt(ix(:)==1);
    dx = diff(sts);
    dx(sign(dx)==-1) =[];
    [isi1(it,:)] = hist(dx,0:501);
    
    ix = x2';
    sts = xt(ix(:)==1);
    dx = diff(sts);
    dx(sign(dx)==-1) =[];
    [isi2(it,:)] = hist(dx,0:501);
    
    if sum(sum(x1,2)>8)~=0
        [s1,spkFx,r1] = mtspectrumpb(x1(sum(x1,2)>8,:)',params1);
        S1(it,:) = s1./r1;
        R1(it) = r1;
    end;
    if sum(sum(x2,2)>8)~=0
        [s2,~,r2] = mtspectrumpb(x2(sum(x2,2)>8,:)',params1);
        S2(it,:) = s2./r2;
        R2(it) = r2;
    end;
                
    ix = [];
    for jt = 1:size(x1,1)
        ix(jt) = jt+size(x1,1)*(jt-1);
    end;

    if sum(sum(x1,2)>8)~=0
        xc1 = xcorr(x1',lag);   
        xc1 = mean(xc1(:,ix(sum(x1,2)>8)),2);
        xc1(501) = min(xc1);
        XC1(it,:) = xc1;
    end;
    
    if sum(sum(x2,2)>8)~=0
        xc2 = xcorr(x2',lag);
        xc2 = mean(xc2(:,ix(sum(x2,2)>8)),2);
        xc2(501) = min(xc2);
        XC2(it,:) = xc2;
    end;
    
    fprintf('\n');
end;

%%
lag = 500;
cnt1 = 0;
sta = {};
for curLFP = 1:length(lfpPoolLabel)
    
    [lfp] = LFPpoolH{curLFP};
    [nsmp,ntrl] = size(lfp);
    
    unitIx = find(strcmp(lfpPoolLabel(curLFP),unitPoolLabel));
    
    for curUnit = 1:length(unitIx)
        cnt1 = cnt1+1;
        spk = spkPool{unitIx(curUnit)};
        if (delIxPool{curLFP}~=0)
            spk(delIxPool{curLFP},:) = [];
        end;
        spk = spk(ixHPool{curLFP},:)';
        
        cnt2 = 0;       
        for curTrl = 1:ntrl
            tsIx = find(spk(:,curTrl)==1);
            for curTs = 1:length(tsIx)
                if (tsIx(curTs)-lag>0) && (tsIx(curTs)+lag <= nsmp)
                    cnt2 = cnt2+1;
                    x = lfp(tsIx(curTs)-lag:tsIx(curTs)+lag,curTrl);
                    x = x -mean(x);
                    sta{cnt1}(:,cnt2) = x;
                end;
            end;
        end;                       
    end;
end;

%%
foi = 2.^((0:53)./8);

pvalR = zeros(length(sta),length(foi));
zval = zeros(length(sta),length(foi));
parfor kt = 1:length(sta)    
    
    [nsmp,ntrl] = size(sta{kt});
    ix = 1:nsmp;
    concat = zeros(1,nsmp*ntrl);
    for jt = 1:ntrl
        concat(ix) = sta{kt}(:,jt);
        ix = ix+nsmp;
    end;
    [phi,~ ] = extractPhaseLFP({concat'},foi,1e3,[1 0]);
    phi{1} = reshape(phi{1},[length(foi) nsmp ntrl]);
    
    p = zeros(1,length(foi));
    z = zeros(1,length(foi));
    for it = 1:length( foi )
        [p(it),z(it)] =  circ_rtest(squeeze(phi{1}(it,501,:)));
    end;
    pvalR(kt,:) = p;
    zval(kt,:) = z;
    
end;

[hcR, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvalR,0.05);

%%
sigIx = [];
for it = 1:size(pvalR,1)
    
    %if any(pvalR(it,:)<0.001/length(foi))
    if any(hcR(it,:)==1)
        sigIx = [sigIx it];
    end;
    
end;

for it = 1:length(sigIx)
    figure;subplot(121);plot(foi,zval(sigIx(it),:));subplot(122);plot(-lag:lag,mean(sta{sigIx(it)},2));
end;

%% 
pow = cell(1,length(sta));
for it = 1:length(sta)
    [~,pow{it}] = extractPhaseLFP({sta{it}},2.^((0:53)./8),1e3,[0 1]);
end;

%%
[staS,staFx] = mtspectrumc(sta',params1);
figure;
subplot(121);
pcolor(-lag:lag,2.^((0:53)./8),10*log10(pow{1})-10*log10(mean(pow{1},2)*ones(1,size(pow{1},2))));
%pcolor(-lag:lag,2.^((0:53)./8),10*log10(pow{1}));
%y = log10(staS);
%x = log10(staFx)';
%b = regress(y,[ones(size(x)) x]);
%yp = [ones(size(x)) x]*b;
%c = 10.^(y-yp);
%plot(staFx,c);
subplot(122);
plot(-lag:lag,mean(sta,1));
        %grid on;