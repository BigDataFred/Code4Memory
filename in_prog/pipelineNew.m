%%
restoredefaultpath;
addpath('~rouxf/prj/Bham/code/mcode/project_EM/');
addpath(genpath('~rouxf/prj/Bham/code/mcode/utils/'));
addpath(genpath('~/tbx/chronux_2_11/spectral_analysis/'));
addpath('~/tbx/fieldtrip-20170618/');
ft_defaults;

%% recruit workers for parallel computing
if isempty(gcp('nocreate'))
    parpool(36,'SpmdEnabled', false)
end;

%%
Fs = 1e3;
[lpf] = fir1(3*fix(Fs/0.5),1/(Fs/0.5),'high');

%%
pID = 'P08';
expMode = 'cnEM';

%%
%p2d = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',pID,'/',expMode,'/'];
p2d = ['/home/rouxf/tmp/resultsC4M/',pID,'/'];

%sesh = dir(p2d);
%sesh(1:2) = [];
%chck = regexp({sesh(:).name},'\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}');
%selIx= []; cnt = 0;for it = 1:length(chck); if (~isempty(chck{it})) && (sesh(it).isdir);cnt = cnt+1; selIx(cnt) = it;end;end;
%sesh = {sesh(selIx).name}'

sesh = dir([p2d,pID,'_',expMode,'_','*_lfpDataStimLockedSegmenteddownsampled.mat']);
chck = regexp({sesh(:).name},'\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}');
dum = cell(1,length(chck));
for it = 1:length( dum )
    dum(it) = {sesh(it).name(chck{it}:chck{it}+18)};
end;
sesh = dum'

%%
[avgSxx1] = cell(length(sesh),1);
[avgSxx2] = cell(length(sesh),1);
[avgSxx3] = cell(length(sesh),1);
[avgSxx4] = cell(length(sesh),1);
[avgFR] = cell(length(sesh),1);
[semFR] = cell(length(sesh),1);

%%
for curSesh = 1:length(sesh)
           
    %%
    LFPfn = [pID,'_',expMode,'_',sesh{curSesh},'_lfpDataStimLockedSegmenteddownsampled.mat'];
    [~,fn,ext] = fileparts(LFPfn);    
    saveName = [fn,'_IEDdetection',ext];
    
    %[lfpDat] = load( [p2d,sesh{curSesh},'/',LFPfn] )
    [lfpDat] = load( [p2d,LFPfn] )
    
    [lfpDat] = backWardComp(lfpDat);
    
    tIx = find(lfpDat.trlTime >=-1 & lfpDat.trlTime <=5);
    [delIx,selIx] = detectIED2trlIx(lfpDat,tIx);
    %save([p2d,sesh{curSesh},'/',saveName],'delIx','selIx');
    save([p2d,saveName],'delIx','selIx');
    
end;

%%
BFlab1 = {};
[LFPpool1] = {};
[LFPpool2] = {};
[LFPpool3] = {};
dt2 = -1e3:250:5e3;
for curSesh = 1%:length(sesh)
    
    %%
    p2SpkFile = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',pID,'/spkRes/'];
    spkFile = [ pID,'_',expMode,'_',sesh{curSesh},'_spkResp.mat' ];
    [spkDat] = load([p2SpkFile,spkFile])   
    
    %%
    selIx = spkDat.sigSel;
    pct = zeros(length(selIx),1);
    for curMW = 1:length(selIx)
        pct(curMW) = sum(spkDat.isiH{selIx(curMW)}(1:4))/sum(spkDat.isiH{selIx(curMW)});
    end;
    mxT = zeros(length(selIx),1);
    for curMW = 1:size(spkDat.fr,1)    
        [~,mIx] = max(spkDat.fr(curMW,:));
        mxT(curMW) = dt2(mIx);
    end;
    spkDat.sigSel = spkDat.sigSel(pct <0.03 & mxT>2e3);
    
    %%
    LFPfn = [pID,'_',expMode,'_',sesh{curSesh},'_lfpDataStimLockedSegmenteddownsampled.mat'];
    [~,fn,ext] = fileparts(LFPfn);
    
    [lfpDat] = load( [p2d,sesh{curSesh},'/',LFPfn] )
    
    %%
    [lfpDat] = backWardComp(lfpDat);
    
    %%    
    IEDfn = [fn,'_IEDdetection',ext];    
    load([p2d,sesh{curSesh},'/',IEDfn]);
    
    for curTrl = 1:length(selIx)
        for curMW = 1:length( selIx{curTrl} )
            lfpDat.LFPseg{selIx{curTrl}(curMW)}(:,delIx{curTrl}) = [];
        end;
    end;
    
    %%
    [lfpDat] = sigProc4TimeFreqAnalysis(lfpDat);
    
    %%
    
    %%
    [hitIdx] = find(ismember(sort([lfpDat.missIdx;lfpDat.hitIdx]),lfpDat.hitIdx));
    [missIdx] = find(ismember(sort([lfpDat.missIdx;lfpDat.hitIdx]),lfpDat.missIdx));
    hIx = cell(1,length(lfpDat.chanLab));
    mIx = cell(1,length(lfpDat.chanLab));
    for it = 1:length(delIx)
        for jt = 1:length(selIx{it})
            tmp1 = hitIdx;
            tmp2 = missIdx;
            tmp1(ismember(hitIdx,delIx{it})) = [];
            tmp2(ismember(missIdx,delIx{it})) = [];
            hIx{selIx{it}(jt)} = find(ismember(sort([tmp1;tmp2]),tmp1));
            mIx{selIx{it}(jt)} = find(ismember(sort([tmp1;tmp2]),tmp2));
        end;
    end;
    
    %%
    [lab,~,~,~] = extractMWLabel(lfpDat.chanLab(spkDat.spkSelIdx(spkDat.sigSel)));
    BFlab1 = [BFlab1 lfpDat.chanLab(spkDat.spkSelIdx(spkDat.sigSel))];
    
    %tIx = find(lfpDat.trlTime >=2 & lfpDat.trlTime <=4);
    
    %[LFPpowDat1] = spectralAnalysisEM2(lfpDat,spkDat.spkSelIdx(spkDat.sigSel),[],5,tIx,true,true,true);
    
    [LFPpowDat2] = spectralAnalysisEM2(lfpDat,spkDat.spkSelIdx(spkDat.sigSel),hIx(spkDat.spkSelIdx(spkDat.sigSel)),5,[],true,true,true);
    %[LFPpowDat3] = spectralAnalysisEM2(lfpDat,spkDat.spkSelIdx(spkDat.sigSel),mIx(spkDat.spkSelIdx(spkDat.sigSel)),5,[],true,true,true);
    
    %[LFPpool1] = [LFPpool1;LFPpowDat1];
    [LFPpool2] = [LFPpool2;LFPpowDat2];
    %[LFPpool3] = [LFPpool3;LFPpowDat3];
    
    %%        
%     [avgSxx1{curSesh},avgSxx2{curSesh},BFid1] = avgLFPPowData(LFPpowDat1,BFlab1);
    
    %%
%     [BFlab2,~,~,~] = extractMWLabel(lfpDat.chanLab);
%     [chanIx] = filterChanels4Label(BFlab2,'Hipp');
%     
%     BFlab2 = BFlab2(chanIx);
%     
%     [statInfo] = baseVSpostDepTtest4LfpPower(lfpDat,[4 12],[2 4],chanIx);
    
%     %%      
%     BFlab2 = BFlab2([statInfo.selIxN;statInfo.selIxP]);
%     
%     tIx = find(lfpDat.trlTime >=2 & lfpDat.trlTime <=4);
%     [LFPpowDat2] = spectralAnalysisEM2(lfpDat,[statInfo.selIxN;statInfo.selIxP],[],5,tIx,true,true);
%     
%     [LFPpool2] = [LFPpool2;LFPpowDat1];    
% 
%     %%        
%     [avgSxx3{curSesh},avgSxx4{curSesh},BFid2] = avgLFPPowData(LFPpowDat2,BFlab2);
    
    %%
%     [avgFR{curSesh},semFR{curSesh},~] = avgFR4BFlabel(spkDat,BFlab1);
     
end;

%%
for it = 1:length(BFlab1)
    x = BFlab1(it);
    x = x{1}(1:end-1);
    BFlab1(it) = {x};
end;
BFid = unique(BFlab1);

%%
avgTFR = zeros(length(BFid), length(LFPpool2{1}.txx1), length(LFPpool2{1}.fxx1) );

for it = 1:length(BFid)
    
    ix = find(strcmp(BFlab1,BFid(it)));
    
    cnt = 0;
    for jt = 1:length( ix )        
        if size(LFPpool2{ix(jt)}.Sxx1,3)>20
            cnt = cnt+1;
            Y = normalize( squeeze(mean(LFPpool2{ix(jt)}.Sxx1,3)) );        
            avgTFR(it,:,:) = squeeze(avgTFR(it,:,:))+Y;        
        end;
    end;
    avgTFR(it,:,:) = avgTFR(it,:,:)./cnt;
    
end;

figure;
for it = 1:size(avgTFR,1)
    subplot(size(avgTFR,1),1,it);
    imagesc(LFPpool2{1}.txx1,LFPpool2{1}.fxx1,squeeze( avgTFR(it,:,:) )' );
    axis xy;xlim([-1 4]);
    title(BFid(it));colormap jet;
end;
%%
avgTFR = zeros(length(BFid), length(LFPpool2{1}.txx2), length(LFPpool2{1}.fxx2) );

for it = 1:length(BFid)
    
    ix = find(strcmp(BFlab1,BFid(it)));
    
    cnt = 0;
    for jt = 1:length( ix )        
        if size(LFPpool2{ix(jt)}.Sxx2,3)>20
            cnt = cnt+1;            
            Y = normalize(squeeze(mean(LFPpool2{ix(jt)}.Sxx2,3)));        
            avgTFR(it,:,:) = squeeze(avgTFR(it,:,:))+Y;        
        end;
    end;
    avgTFR(it,:,:) = avgTFR(it,:,:)./cnt;
    
end;

figure;
for it = 1:size(avgTFR,1)
    subplot(size(avgTFR,1),1,it);
    imagesc(LFPpool2{1}.txx2,LFPpool2{1}.fxx2,squeeze( avgTFR(it,:,:) )' );
    axis xy;xlim([-1 4]);
    title(BFid(it));colormap jet;
end;

%%
for it = 1:length(LFPpool2)
    if ~isempty(LFPpool2{it}) && ~isempty(LFPpool3{it})
        if ( isfield(LFPpool2{it},'Sxx1') ) && ( isfield(LFPpool3{it},'Sxx1') ) %&& ( size(LFPpool2{it}.Sxx1,3)>20 ) && ( size(LFPpool3{it}.Sxx1,3)>20 )
            
            cnt = cnt+1;
            x1 = squeeze(mean(LFPpool2{it}.Sxx1,3));
            x2 = squeeze(mean(LFPpool3{it}.Sxx1,3));
            Y = normalize( [x1;x2] );
            
            avgTFR1 = avgTFR1 + Y(1:size(x1,1),:);
            avgTFR2 = avgTFR2 + Y(size(x1,1)+1:size(x1,1)+size(x2,1),:);
            
            x1 = squeeze(mean(LFPpool2{it}.Sxx2,3));
            x2 = squeeze(mean(LFPpool3{it}.Sxx2,3));
            Y = normalize( [x1;x2] );
            
            avgTFR3 = avgTFR3 + Y(1:size(x1,1),:);
            avgTFR4 = avgTFR4 + Y(size(x1,1)+1:size(x1,1)+size(x2,1),:);
            
        end;
    end;
end;

%%
avgTFR1 = zeros(size(LFPpool2{1}.Sxx1,1),size(LFPpool2{1}.Sxx1,2));
avgTFR2 = zeros(size(LFPpool3{1}.Sxx1,1),size(LFPpool3{1}.Sxx1,2));
avgTFR3 = zeros(size(LFPpool2{1}.Sxx2,1),size(LFPpool2{1}.Sxx2,2));
avgTFR4 = zeros(size(LFPpool3{1}.Sxx2,1),size(LFPpool3{1}.Sxx2,2));

figure;
cnt = 0;
for it = 1:length(LFPpool2)
    if ~isempty(LFPpool2{it}) && ~isempty(LFPpool3{it})
        if ( isfield(LFPpool2{it},'Sxx1') ) && ( isfield(LFPpool3{it},'Sxx1') ) %&& ( size(LFPpool2{it}.Sxx1,3)>20 ) && ( size(LFPpool3{it}.Sxx1,3)>20 )
            cnt = cnt+1;
            x1 = squeeze(mean(LFPpool2{it}.Sxx1,3));
            x2 = squeeze(mean(LFPpool3{it}.Sxx1,3));
            Y = normalize( [x1;x2] );
            
            avgTFR1 = avgTFR1 + Y(1:size(x1,1),:);
            avgTFR2 = avgTFR2 + Y(size(x1,1)+1:size(x1,1)+size(x2,1),:);
            
            x1 = squeeze(mean(LFPpool2{it}.Sxx2,3));
            x2 = squeeze(mean(LFPpool3{it}.Sxx2,3));
            Y = normalize( [x1;x2] );
            
            avgTFR3 = avgTFR3 + Y(1:size(x1,1),:);
            avgTFR4 = avgTFR4 + Y(size(x1,1)+1:size(x1,1)+size(x2,1),:);
        end;
    end;
end;
avgTFR1 = avgTFR1./cnt;
avgTFR2 = avgTFR2./cnt;
avgTFR3 = avgTFR3./cnt;
avgTFR4 = avgTFR4./cnt;

figure;
subplot(4,2,[1 3 5]);
imagesc(LFPpool2{1}.txx2,LFPpool2{1}.fxx2,avgTFR3');
axis xy;ca = caxis;xlim([-1 5]);
colorbar;
subplot(4,2,[2 4 6]);
imagesc(LFPpool2{1}.txx2,LFPpool2{1}.fxx2,avgTFR4');
axis xy;caxis(ca);xlim([-1 5]);
colorbar;
subplot(4,2,7);
imagesc(LFPpool2{1}.txx1,LFPpool2{1}.fxx1,avgTFR1');
axis xy;ca = caxis;xlim([-1 5]);
colorbar;
subplot(4,2,8);
imagesc(LFPpool2{1}.txx1,LFPpool2{1}.fxx1,avgTFR2');
axis xy;caxis(ca);xlim([-1 5]);
colorbar;

%%
tIx = find(LFPpool1{1}.txx1 >=-2 & LFPpool1{1}.txx1 <=0);
for curMW = 1:length(LFPpool1)
    
    x1 =  LFPpool1{curMW}.Sxx1 ;    
    x2 = mean( LFPpool1{curMW}.Sxx1(tIx,:,:),1);
    [k,p,~] = size(x1);
    empT = zeros(k,p);
        
    for it = 1:k
        for jt = 1:p
            
            d = x1(it,jt,:) - x2(1,jt,:);
                        
            empT(it,jt) = mean(d)/(std(d)/length(d));
        end;
    end;
    figure;
    imagesc(LFPpool2{1}.txx1,LFPpool2{1}.fxx1,empT');
    axis xy;xlim([0 5]);
    
end;

%%
d = [];
cnt = 0;
for curMW = 1:length(LFPpool2)
    
    if ~isempty(LFPpool2{curMW}) && ~isempty(LFPpool3{curMW})
        if size(LFPpool2{curMW}.Sxx1,3) >20 && size(LFPpool3{curMW}.Sxx1,3)>20
            cnt = cnt+1;
            x1 =  squeeze(mean(LFPpool2{curMW}.Sxx1,3)) ; %Hits
            x2 =  squeeze(mean(LFPpool3{curMW}.Sxx1,3)) ; %Misses       
            y = normalize([x1;x2]);
            d(cnt,:,:) = y(1:size(x1,1),:)-y(size(x1,1)+1:size(x1,1)+size(x2,1),:);
        end;
    end;
end;


[k,p,~] = size(d);
empT = zeros(k,p);
for it = 1:k
    for jt = 1:p
        
        [~,~,~,stat] = ttest(squeeze(d(:,it,jt)),0);
        empT(it,jt) = stat.tstat;
        
    end;
end;

%%
[X1] = zeros(size(LFPpool2,1),size(LFPpool2{1}.Sxx1,1),size(LFPpool2{1}.Sxx1,2));
[X2] = zeros(size(LFPpool3,1),size(LFPpool3{1}.Sxx1,1),size(LFPpool3{1}.Sxx1,2));

cnt = 0;
for it = 1:size( LFPpool2,1 )
    
    if size(LFPpool2{it}.Sxx1,3) >20 && size(LFPpool3{it}.Sxx1,3)>20
        cnt = cnt+1;
        X1(cnt,:,:) = squeeze(mean(LFPpool2{it}.Sxx1,3));% Hits 
        X2(cnt,:,:) = squeeze(mean(LFPpool3{it}.Sxx1,3));% Misses
    end;
end;

[n,p,~] = size(LFPpool2{1}.Sxx1);
empT = zeros(n,p);
for it = 1:n
    for jt = 1:p
        
        d = X1(:,it,jt)-X2(:,it,jt); % Hits - Misses
        df = sqrt(length(d));
        empT(it,jt) = mean(d)/(std(d)/df);
    end;
end;
figure;
imagesc(LFPpool2{1}.txx1,LFPpool2{1}.fxx1,empT)
axis xy;

%%
tIx1 = find(LFPpool1{1}.txx1 >= -2 &  LFPpool1{1}.txx1 <=0);   
tIx2 = find(LFPpool1{1}.txx1 >=  2 &  LFPpool1{1}.txx1 <=4);
fIx = find(LFPpool1{1}.fxx1 >= 4 &  LFPpool1{1}.fxx1 <= 11);

[avgPow1] = zeros( size(LFPpool1{1}.Sxx1) );
[avgPow2] = zeros( size(LFPpool1{1}.Sxx2) );
[pow] = zeros( length( LFPpool1 ),2 );
for curMW = 1:length( LFPpool1 )
        
    Y1 = LFPpool2{curMW}.Sxx1;%normalize( LFPpool1{curMW}.Sxx1 );
    Y2 = LFPpool2{curMW}.Sxx2;%normalize( LFPpool1{curMW}.Sxx2 );
    
    x1 = mean(mean(Y1(tIx1,fIx)));
    x2 = mean(mean(Y2(tIx2,fIx)));
    
    pow(curMW,:) = [x1 x2];
        
    [avgPow1] = avgPow1+Y1;        
    [avgPow2] = avgPow2+Y2;
end;
avgPow1 = avgPow1./length( LFPpool1 );
avgPow2 = avgPow2./length( LFPpool1 );

figure;
subplot(221);
imagesc(LFPpool1{1}.txx2,LFPpool1{1}.fxx2,avgPow2');
axis xy;colormap jet;
xlim([-2 5]);
subplot(222);
subplot(223);
imagesc(LFPpool1{1}.txx1,LFPpool1{1}.fxx1,avgPow1');
axis xy;colormap jet;
xlim([-2 5]);
subplot(224);
hold on;
plot(ones(1,size(pow,1)),pow(:,1),'bo');
plot(2*ones(1,size(pow,1)),pow(:,2),'ro');
for it = 1:size(pow,1)
    plot([1 2],[pow(it,1) pow(it,2)],'-','Color',[.75 .75 .75]);
end;
xlim([0 3]);

%%
tIx1 = find(LFPpool1{1}.txx1 >= -2 &  LFPpool1{1}.txx1 <=0);   
tIx2 = find(LFPpool1{1}.txx1 >=  2 &  LFPpool1{1}.txx1 <=4);
fIx = find(LFPpool1{1}.fxx1 >= 4 &  LFPpool1{1}.fxx1 <= 11);

[avgPow1] = zeros( size(LFPpool1{1}.Sxx1) );
[avgPow2] = zeros( size(LFPpool1{1}.Sxx2) );
[pow1] = zeros( length(LFPpool1),2 );
for curMW = 1:length( LFPpool1 )
    
    
    Y1 = LFPpool3{curMW}.Sxx1;%(LFPpool1{curMW}.Sxx1-min(min([LFPpool1{curMW}.Sxx1;LFPpool2{curMW}.Sxx1])))./(max(max([LFPpool1{curMW}.Sxx1;LFPpool2{curMW}.Sxx1]))-min(min([LFPpool1{curMW}.Sxx1;LFPpool2{curMW}.Sxx1])));
    Y2 = LFPpool3{curMW}.Sxx2;%(LFPpool1{curMW}.Sxx2-min(min([LFPpool1{curMW}.Sxx2;LFPpool2{curMW}.Sxx2])))./(max(max([LFPpool1{curMW}.Sxx2;LFPpool2{curMW}.Sxx2]))-min(min([LFPpool1{curMW}.Sxx2;LFPpool2{curMW}.Sxx2])));
    
    x1 = mean(mean( Y1( tIx1,fIx) ) );
    x2 = mean(mean( Y1( tIx2,fIx) ) );
    
    pow1(curMW,:) = [x1 x2];
    
    [avgPow1] = avgPow1+Y1;        
    [avgPow2] = avgPow2+Y2;
end;
avgPow1 = avgPow1./length( LFPpool1 );
avgPow2 = avgPow2./length( LFPpool1 );

figure;
subplot(221);
imagesc(LFPpool1{1}.txx2,LFPpool1{1}.fxx2,avgPow2');
axis xy;colormap jet;
xlim([-2 5]);
subplot(222);
subplot(223);
imagesc(LFPpool1{1}.txx1,LFPpool1{1}.fxx1,avgPow1');
axis xy;colormap jet;
xlim([-2 5]);
subplot(224);
hold on;
plot(ones(1,size(pow1,1)),pow1(:,1),'bo');
plot(2*ones(1,size(pow1,1)),pow1(:,2),'ro');
for it = 1:size(pow1,1)
    plot([1 2],[pow1(it,1) pow1(it,2)],'-','Color',[.75 .75 .75]);
end;
xlim([0 3]);

%%
tIx1 = find(LFPpool2{1}.txx1 >= -2 &  LFPpool2{1}.txx1 <=0);   
tIx2 = find(LFPpool2{1}.txx1 >=  2 &  LFPpool2{1}.txx1 <=4);
fIx = find(LFPpool2{1}.fxx1 >= 4 &  LFPpool2{1}.fxx1 <= 11);

[avgPow1] = zeros( size(LFPpool2{1}.Sxx1) );
[avgPow2] = zeros( size(LFPpool2{1}.Sxx2) );
[pow1] = zeros( length(LFPpool2),2 );
for curMW = 1:length( LFPpool2 )
    
    Y1 = (LFPpool2{curMW}.Sxx1-min(min([LFPpool1{curMW}.Sxx1;LFPpool2{curMW}.Sxx1])))./(max(max([LFPpool1{curMW}.Sxx1;LFPpool2{curMW}.Sxx1]))-min(min([LFPpool1{curMW}.Sxx1;LFPpool2{curMW}.Sxx1])));
    Y2 = (LFPpool2{curMW}.Sxx2-min(min([LFPpool1{curMW}.Sxx2;LFPpool2{curMW}.Sxx2])))./(max(max([LFPpool1{curMW}.Sxx2;LFPpool2{curMW}.Sxx2]))-min(min([LFPpool1{curMW}.Sxx2;LFPpool2{curMW}.Sxx2])));
    
    x1 = mean(mean( Y1( tIx1,fIx) ) );
    x2 = mean(mean( Y1( tIx2,fIx) ) );
    
    pow1(curMW,:) = [x1 x2];
    
    [avgPow1] = avgPow1+Y1;
    
    Y2 = LFPpool2{curMW}.Sxx2;
    Y2 = normalize(Y2);
    [avgPow2] = avgPow2+Y2;
end;
avgPow1 = avgPow1./length( LFPpool2 );
avgPow2 = avgPow2./length( LFPpool2 );

figure;
subplot(221);
imagesc(LFPpool2{1}.txx2,LFPpool2{1}.fxx2,avgPow2');
axis xy;colormap jet;
xlim([-2 5]);
subplot(222);
subplot(223);
imagesc(LFPpool2{1}.txx1,LFPpool2{1}.fxx1,avgPow1');
axis xy;colormap jet;
xlim([-2 5]);
subplot(224);
hold on;
plot(ones(1,size(pow1,1)),pow1(:,1),'bo');
plot(2*ones(1,size(pow1,1)),pow1(:,2),'ro');
for it = 1:size(pow1,1)
    plot([1 2],[pow1(it,1) pow1(it,2)],'-','Color',[.75 .75 .75]);
end;
xlim([0 3]);

%%
tIx1 = find(LFPpool1{1}.txx1 >= 2 & LFPpool1{1}.txx1 <= 4);
tIx2 = find(LFPpool1{1}.txx2 >= 2 & LFPpool1{1}.txx2 <= 4);

fIx1 = find(LFPpool1{1}.fxx1 >= 4 &  LFPpool1{1}.fxx1 <= 11);
fIx2 = find(LFPpool1{1}.fxx2 >= 30 &  LFPpool1{1}.fxx2 <= 60);

[avgPow1] = zeros( size(LFPpool1{1}.Sxx1) );
[avgPow2] = zeros( size(LFPpool1{1}.Sxx2) );
[lowPow] = zeros( length(LFPpool1),2 );
[highPow] = zeros( length(LFPpool1),2 );

for curMW = 1:length( LFPpool1 )
        
    x1 = (LFPpool1{curMW}.Sxx1);%Hits
    x2 = (LFPpool1{curMW}.Sxx2);%Hits
    x3 = (LFPpool2{curMW}.Sxx1);%Misses
    x4 = (LFPpool2{curMW}.Sxx2);%Misses
    
    Y1 = normalize(x1-x3);
    Y2 = normalize(x2-x4);   
    
    lowPow(curMW,:) = [mean(mean(x1(tIx1,fIx1),1)) mean(mean(x3(tIx1,fIx1),1))];
    highPow(curMW,:) = [mean(mean(x2(tIx2,fIx2),1)) mean(mean(x4(tIx2,fIx2),1))];
    
    [avgPow1] = avgPow1+Y1;        
    [avgPow2] = avgPow2+Y2;
end;
avgPow1 = avgPow1./length( LFPpool1 );
avgPow2 = avgPow2./length( LFPpool1 );

figure;
subplot(221);
imagesc(LFPpool1{1}.txx2,LFPpool1{1}.fxx2,avgPow2');
axis xy;colormap jet;
subplot(223);
imagesc(LFPpool1{1}.txx1,LFPpool1{1}.fxx1,avgPow1');
axis xy;colormap jet;
subplot(222);
hold on;
plot(ones(1,size(highPow,1)),highPow(:,1),'bo');
plot(2*ones(1,size(highPow,1)),highPow(:,2),'ro');
for it = 1:size(highPow,1)
    plot([1 2],[highPow(it,1) highPow(it,2)],'-','Color',[.75 .75 .75]);
end;
xlim([0 3]);
subplot(224);
hold on;
plot(ones(1,size(lowPow,1)),lowPow(:,1),'bo');
plot(2*ones(1,size(lowPow,1)),lowPow(:,2),'ro');
for it = 1:size(lowPow,1)
    plot([1 2],[lowPow(it,1) lowPow(it,2)],'-','Color',[.75 .75 .75]);
end;
xlim([0 3]);

%%
for curSesh = 1:length( sesh )
    
    figure;
    n = length(avgFR{curSesh});
    for it = 1:n
        subplot(1,n,it);
        jbfill(-1e3:250:5e3,avgFR{curSesh}{it}-semFR{curSesh}{it},avgFR{curSesh}{it}+semFR{curSesh}{it},[.9 0 0],[.9 0 0],0,.5);
        hold on;
        plot(-1e3:250:5e3,avgFR{curSesh}{it},'r-','LineWidth',3);
        plot([-1e3 5e3],[0 0],'k--');
        %plot([0 0],[min(avgFR{curSesh}{it}-semFR{curSesh}{it}) max(avgFR{curSesh}{it}+semFR{curSesh}{it})],'k');
        %plot([2e3 2e3],[min(avgFR{curSesh}{it}-semFR{curSesh}{it}) max(avgFR{curSesh}{it}+semFR{curSesh}{it})],'k');
        %title(BFid1{it});
        axis tight;xlim([-1e3 5e3]);        
    end;
    
    %%
    figure;
    hold on;
    for it = 1:length(LFPpowDat1)
        subplot(3,ceil(length(LFPpowDat1)/3),it);
        hold on;
        imagesc((LFPpowDat1{it}.txx2).*1e3,LFPpowDat1{it}.fxx2,(LFPpowDat1{it}.Sxx2)');axis xy;colormap jet;
        plot([0 0],[min(LFPpowDat1{it}.fxx2) max(LFPpowDat1{it}.fxx2)],'w-');
        plot([2 2].*1e3,[min(LFPpowDat1{it}.fxx2) max(LFPpowDat1{it}.fxx2)],'w-');
        axis tight;xlim([-1.5 5].*1e3);
        %title(lfpDat.chanLab{spkDat.spkSelIdx(it)});
    end;
    
    figure;
    n = length(avgSxx1{curSesh});
    for it = 1:n
        subplot(1,n,it);
        hold on;
        imagesc((LFPpowDat1{1}.txx2).*1e3,LFPpowDat1{1}.fxx2,avgSxx2{curSesh}{it}');axis xy;colormap jet;
        plot([0 0],[min(LFPpowDat1{1}.fxx2) max(LFPpowDat1{1}.fxx2)],'w-');
        plot([2 2].*1e3,[min(LFPpowDat1{1}.fxx2) max(LFPpowDat1{1}.fxx2)],'w-');
        axis tight;xlim([-1.5 5].*1e3);
        %title(BFid1{it});
        colorbar;
    end;
    
    %%
    figure;
    hold on;
    for it = 1:length(LFPpowDat1)
        subplot(3,ceil(length(LFPpowDat1)/3),it);
        hold on;
        imagesc((LFPpowDat1{it}.txx1).*1e3,LFPpowDat1{it}.fxx1,(LFPpowDat1{it}.Sxx1)');axis xy;colormap jet;
        plot([0 0],[min(LFPpowDat1{it}.fxx1) max(LFPpowDat1{it}.fxx1)],'w-');
        plot([2 2].*1e3,[min(LFPpowDat1{it}.fxx1) max(LFPpowDat1{it}.fxx1)],'w-');
        axis tight;xlim([-1.5 5].*1e3);
        %title(lfpDat.chanLab{spkDat.spkSelIdx(it)});
    end;
    
    figure;
    n = length(avgSxx1{curSesh});
    for it = 1:n
        subplot(1,n,it);
        hold on;
        imagesc((LFPpowDat1{1}.txx1).*1e3,LFPpowDat1{1}.fxx1,avgSxx1{curSesh}{it}');axis xy;colormap jet;
        plot([0 0],[min(LFPpowDat1{1}.fxx1) max(LFPpowDat1{1}.fxx1)],'w-');
        plot([2 2].*1e3,[min(LFPpowDat1{1}.fxx1) max(LFPpowDat1{1}.fxx1)],'w-');
        axis tight;xlim([-1.5 5].*1e3);
        %title(BFid1{it});
        colorbar;
    end;
    
    %%
    chanIx = [statInfo.selIxN;statInfo.selIxP];
    figure;
    hold on;
    for it = 1:length(LFPpowDat2)
        subplot(3,ceil(length(LFPpowDat2)/3),it);
        hold on;
        imagesc((LFPpowDat2{it}.txx2).*1e3,LFPpowDat2{it}.fxx2,(LFPpowDat2{it}.Sxx2)');axis xy;colormap jet;
        plot([0 0],[min(LFPpowDat2{it}.fxx2) max(LFPpowDat2{it}.fxx2)],'w-');
        plot([2 2].*1e3,[min(LFPpowDat2{it}.fxx2) max(LFPpowDat2{it}.fxx2)],'w-');
        axis tight;xlim([-1.5 5].*1e3);
        %title(lfpDat.chanLab{chanIx(it)});
    end;
    
    figure;
    n = length(avgSxx3{curSesh});
    for it = 1:n
        subplot(1,n,it)
        hold on;
        imagesc((LFPpowDat2{1}.txx2).*1e3,LFPpowDat2{1}.fxx2,avgSxx4{curSesh}{it}');axis xy;colormap jet;
        plot([0 0],[min(LFPpowDat2{1}.fxx2) max(LFPpowDat2{1}.fxx2)],'w-');
        plot([2 2].*1e3,[min(LFPpowDat2{1}.fxx2) max(LFPpowDat2{1}.fxx2)],'w-');
        axis tight;xlim([-1.5 5].*1e3);
        %title(BFid2{it});
        colorbar;
    end;
    
    %%
    figure;
    hold on;
    for it = 1:length(LFPpowDat2)
        subplot(3,ceil(length(LFPpowDat2)/3),it);
        hold on;
        imagesc((LFPpowDat2{it}.txx1).*1e3,LFPpowDat2{it}.fxx1,(LFPpowDat2{it}.Sxx1)');axis xy;colormap jet;
        plot([0 0],[min(LFPpowDat2{it}.fxx1) max(LFPpowDat2{it}.fxx1)],'w-');
        plot([2 2].*1e3,[min(LFPpowDat2{it}.fxx1) max(LFPpowDat2{it}.fxx1)],'w-');
        axis tight;xlim([-1.5 5].*1e3);
        %title(lfpDat.chanLab{chanIx(it)});
    end;
    
    figure;
    n = length(avgSxx4{curSesh});
    for it = 1:n
        subplot(1,n,it)
        hold on;
        imagesc((LFPpowDat2{1}.txx1).*1e3,LFPpowDat2{1}.fxx1,avgSxx3{curSesh}{it}');axis xy;colormap jet;
        plot([0 0],[min(LFPpowDat2{1}.fxx1) max(LFPpowDat2{1}.fxx1)],'w-');
        plot([2 2].*1e3,[min(LFPpowDat2{1}.fxx1) max(LFPpowDat2{1}.fxx1)],'w-');
        axis tight;xlim([-1.5 5].*1e3);
        %title(BFid2{it});
        colorbar;
    end;
    
end;