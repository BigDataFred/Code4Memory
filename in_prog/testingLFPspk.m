[pID] = 'P09';
[expMode] = 'fVSpEM';

%%
p2dOrig = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',pID,'/',expMode,'/'];
p2d     = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',pID,'/res/'];
savePath = p2d;

%%
sesh = dir(p2dOrig);
[sesh] = sesh([sesh.isdir]==1);
chck = regexp({sesh.name},'\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}');
cnt = 0;
sel = [];
for it = 1:length( chck )
    if ~isempty( chck{it} )
        cnt = cnt+1;
        sel(cnt) = it;
    end;
end;
[sesh] = {sesh(sel).name}'

%%
for seshIt =1%:3
    
    files= dir([p2d,pID,'_',expMode,'_',sesh{seshIt},'_*_comodulogramDat.mat']);
    
    load([p2d,pID,'_',expMode,'_',sesh{seshIt},'_preprocLFP.mat']);
    
    p2d2 = [p2dOrig,sesh{seshIt},'/'];
    spkDat = dir([p2d2,pID,'_',expMode,'_',sesh{seshIt},'_spkDataStimLockedSegmented.mat']);
    load([p2d2,spkDat.name])
    
    %%
    MWlab = cell(1,length(files));
    BFlab = cell(1,length(files));
    for it = 1:length( BFlab )
        
        startIndex = regexp(files(it).name,'\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}_')+20;
        stopIndex = findstr(files(it).name,'comodulogram')-2;
        d = files(it).name(startIndex:stopIndex);
        MWlab{it} = d;
        
        d = d(1:regexp(d,'\d{1,2}')-1);
        
        BFlab{it} = d;
        
    end;
    BFid = unique(BFlab);
    
    hemLab = cell(1,length(BFid));
    for it  = 1:length( BFid )
        hemLab{it} = BFid{it}(end);
    end;
    [~,sIdx] = sort(hemLab);
    BFid = BFid(sIdx);
    
    %%
    Fs = 1e3;
    dt = -1e3:4e3;
    dx = 500;
    
    TW = 2;
    k = floor(2*TW-1);
    
    params                  = [];
    params.Fs               =Fs;
    params.pad              = 4;
    params.fpass            = [0 30];
    params.tapers           = [TW k];
    params.trialave         = 1;
    
    lpf = fir1(floor(3*(Fs/30)),30/(Fs/2),'low');
    
    
    uC = {'r','b','k','g'};
    
    MI = {};
    S1 = {};
    S2 = {};
    S3 = {};
    S4 = {};
    %PTA = {};
    sigIdx = [];
    cnt1 = 0;
    cnt2 = 0;
    for it = 1:length(MWlab)
        
        fprintf([num2str(it),'/',num2str(length(MWlab))]);
        
        [selIx] = find(strcmp(lfpDat.chanLab,MWlab{it}));
        [ lfx ] = lfpDat.LFPseg{selIx}(lfpDat.trlTime>=-1 & lfpDat.trlTime<=4,:);
        
        [ lfx ] = filtfilt(lpf,1,lfx);% low pass filter LFP <30Hz
        
        fn = dir([p2d,pID,'_',expMode,'_',sesh{seshIt},'_',MWlab{it},'_comodulogramDat.mat']);
        [dat] = load([p2d,fn.name]);
        
        file2= dir([p2d,pID,'_',expMode,'_',sesh{seshIt},'_',MWlab{it},'_spectralData.mat']);
        [dat2] = load([p2d,file2.name]);
        
        % if there is a sig theta response
        if (strcmp(lfpDat.chanLab(selIx),MWlab{it})) && (dat2.pvalPowL<dat.alphaTrsh) %&& any(any(dat.pval<dat.alphaTrsh))
            cnt1 = cnt1+1;
            sigIdx(cnt1) = it;
            
            S1{cnt1} = dat2.S1;
            f1 = dat2.f1;
            S2{cnt1} = dat2.S2;
            f2 = dat2.f2;
            
            S3{cnt1} = dat2.S3;
            t3 = dat2.t3;
            f3 = dat2.f3;
            S4{cnt1} = dat2.S4;
            t4 = dat2.t4;
            f4 = dat2.f4;
            
            MI{cnt1} = dat.MI.*(dat.pval<dat.alphaTrsh);
            afoi = dat.afoi;
            pfoi = dat.pfoi;
            
            %for jt = 1:length(sortedSpikes)
            [spkDat] = sortedSpikes{selIx};%{selIx};% extract the spikeSorted data for the MW of interest
            [cID]    = unique(spkDat.assignedClusterSeg);% cluster labels asigned to each unit
            [trlID]  = lfpDat.hitIdx(selTrl{selIx});% trial-labels for each spike time
            
            if (~isempty(cID))
                for kt = 1:length(cID)
                    
                    [trl] =spkDat.trl(spkDat.assignedClusterSeg == cID(kt));
                    [ts] = spkDat.SpikeTimesSeg(spkDat.assignedClusterSeg == cID(kt)).*1e3;
                    trl = trl(ts>=-1e3 & ts <=4e3);% truncate trial-labels
                    ts  = ts(ts>=-1e3 & ts <=4e3);% truncate spike-times
                    
                    pta = [];
                    cnt3 = 0;
                    for lt = 1:length(trlID)
                        tx = ts(trl == trlID(lt));
                        [n,~] = hist(tx,dt);
                        spkIx = find(n~=0);
                        if (~isempty(spkIx))
                            for mt = 1:length(spkIx)
                                if ( (spkIx(mt)-dx>0) && (spkIx(mt)+dx < size(lfx,1)) )
                                    cnt3 = cnt3+1;
                                    d = lfx( spkIx(mt)-dx:spkIx(mt)+dx, lt );
                                    %PTA{cnt1}(cnt3,:) = d;%
                                    pta(cnt3,:) = d;
                                end;
                            end;
                        end;
                    end;
                    
                    if (size(pta,1)>50)
                        h1 = figure;
                        subplot(4,4,[1:2 5:6]);
                        imagesc(-dx:dx,1:size(pta,1),pta);
                        axis xy;
                        title(MWlab{it});
                        
                        subplot(4,4,[9:10 13:14]);
                        hold on;
                        M = mean(pta,1);
                        SE = std(pta,0,1)/sqrt(size(pta,1)-1);
                        jbfill(-dx:dx,M-SE,M+SE,[.75 .75 .75],[.75 .75 .75],1,.5);
                        hold on;
                        plot(-dx:dx,M,uC{3},'LineWidth',3);
                        plot([-dx dx],[0 0],'r--');
                        axis tight;
                        
                        subplot(4,4,[3:4 7:8]);
                        isi = diff(ts);
                        isi(sign(isi)==-1) = [];
                        [n,~] = hist(isi,[0:1:601]);
                        n(end) = [];
                        plot(0:600,n);
                        xlim([-10 600]);
                        
                        subplot(4,4,[11:12 15:16]);
                        [s1]   = mtspectrumc( (gradient(pta))', params );
                        [s2,f] = mtspectrumc( (gradient(mean(pta,1)))', params );
                        plot(f,s2./s1);
                        
                        %[s,t,f] = mtspecgramc( gradient(pta)', movingwin, params );
                        
                    end;
                    
                end;
            %end;
            end;
        end;
        
        saveName = [pID,'_',expMode,'_',sesh{seshIt},'_',MWlab{it},'_filteredSpectralData.mat'];
        save([savePath,saveName],'S1','f1','S2','f2','S3','t3','f3','S4','t4','f4','MI','afoi','pfoi');
        
        fprintf('\n');
    end;
    
    BFlab = BFlab(sigIdx);
    
    
    %%
    figure;
    a = [];
    b = [];
    ca = [];
    c = 0;
    for it = 1:length( BFid )
        
        [selIdx] = find(strcmp(BFlab,BFid{it}));
        
        AVG = zeros(length(selIdx),size(MI{1},1),size(MI{1},2));
        AVG2 = zeros(length(selIdx),size(S1{1},1));
        AVG3 = zeros(length(selIdx),size(S2{1},1));
        AVG4= zeros(length(selIdx),size(S3{1},1),size(S3{1},2));
        AVG5= zeros(length(selIdx),size(S4{1},1),size(S4{1},2));
        for jt = 1:length(selIdx)
            
            AVG(jt,:,:) = MI{selIdx(jt)};
            AVG2(jt,:) = S1{selIdx(jt)};
            AVG3(jt,:) = S2{selIdx(jt)};
            AVG4(jt,:,:) = S3{selIdx(jt)};
            AVG5(jt,:,:) = S4{selIdx(jt)};
        end;
        
        c = c+1;
        subplot(length(BFid),4,c);
        a = [a gca];
        %SE = std(AVG2,0,1)./sqrt(size(AVG2,1)-1);
        %M = mean(AVG2,1);
        %jbfill(dat2.f1,M-SE,M+SE,[0 0 .9],[0 0 .9],1,.5);
        hold on;
        plot(dat2.f1,AVG2,'k');
        plot(dat2.f2,AVG3,'k');
        plot(dat2.f1,mean(AVG2,1),'r');
        plot(dat2.f2,mean(AVG3,1),'r');
        axis tight;xlim([0 100]);
        title([BFid{it},' n=',num2str(length(selIdx))]);
        
        c = c+1;
        subplot(length(BFid),4,c);
        hold on;
        imagesc(dat2.t3.*1e3,dat2.f3,squeeze(mean( 20*log10(AVG4),1))');
        plot([0 0],[min(dat2.f3) max(dat2.f3)],'w');
        plot([2e3 2e3],[min(dat2.f3) max(dat2.f3)],'w');
        axis xy;axis tight;xlim([-.5 4].*1e3);
        cb = colorbar;
        
        c = c+1;
        subplot(length(BFid),4,c);
        b = [b gca];
        %pcolor(dat.pfoi,dat.afoi, squeeze(mean( AVG ,1)));
        %ca(it,:) = caxis;shading interp;  ylim([40 100]);
        hold on;
        imagesc(dat2.t4.*1e3,dat2.f4,squeeze(mean( 20*log10(AVG5),1))');
        plot([0 0],[min(dat2.f4) max(dat2.f4)],'w');
        plot([2e3 2e3],[min(dat2.f4) max(dat2.f4)],'w');
        axis xy;axis tight;xlim([-.5 4].*1e3);
        cb = colorbar;
        
        c = c+1;
        subplot(length(BFid),4,c);
        b = [b gca];
        pcolor(dat.pfoi,dat.afoi, squeeze(mean( AVG ,1)));
        ca(it,:) = caxis;shading interp;
        cb = colorbar;
        
        colormap jet;
        
    end;
    for it = 1:length(b)
        %caxis(b(it),[min(min(ca)) max(max(ca))]);
    end;
    
end;