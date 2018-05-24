
%%
pID = {'P08'}; %%{'P04','P07','P08','P09'};
expMode = 'fVSpEM';
nrand = 2000;

dt = [-1e3:5e3];
dt2 = [-1.25e3:250:5.25e3];

%%
if isempty(gcp('nocreate'))
    parpool(36,'SpmdEnabled',false);
end;

%%
for pt = 1:length(pID)
    
    p2dOrig = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',pID{pt},'/',expMode,'/'];
    
    files = dir(p2dOrig);
    files(1:2) = [];
    sesh = cell(1,length(files));
    for it = 1:length(files)
        starIx = regexp(files(it).name,'\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}');
        stopIx = starIx+18;
        sesh(it) = {files(it).name(starIx:stopIx)};
    end;
    seshID = unique(sesh);
    
    for seshIt = 1:length(seshID)
        
        [p2d2] = [p2dOrig,seshID{seshIt},'/'];
        [spkDat] = load([p2d2,pID{pt},'_',expMode,'_',seshID{seshIt},'_spkDataStimLockedSegmented.mat'])
        [lfpDat] = load([p2d2,pID{pt},'_',expMode,'_',seshID{seshIt},'_lfpDataStimLockedSegmenteddownsampled.mat'])
        
        try
            load([p2d2,pID{pt},'_',seshID{seshIt},'_LogFile_EMtask_LogDat.mat'])
            logDat
        catch
            p2Log = [p2d2,'log_dat/'];
            logF = dir([p2Log,'*_LogFile_EMtask_LogDat.mat']);
            logDat = load([p2Log,logF.name])
        end;
        
        %%
        [selIx] = sort([lfpDat.missIdx;lfpDat.hitIdx]);
        [x] = logDat.LogDat1.log(selIx,3:4);
        stimCat = cell(length(x),1);
        for it = 1:size(x,1)
            
            stimCat(it) = {[x{it,1}(1) x{it,2}(1)]};
            
        end;
        stimID = unique(stimCat);
        
        %%
        if ~isfield(lfpDat,'dsTrlTime')
            lfpDat.dsTrlTime = lfpDat.trlTime;
        end;
        [ lfp ] = lfpDat.LFPseg{1}(:,find(ismember(lfpDat.trlSel,lfpDat.trlENC)));
        [ lfp ] = lfp(lfpDat.dsTrlTime>=dt(1).*1e-3 & lfpDat.dsTrlTime<=dt(length(dt)).*1e-3,:);
        [lfpDat.dsTrlTime] = lfpDat.dsTrlTime(lfpDat.dsTrlTime>=dt(1).*1e-3 & lfpDat.dsTrlTime<=dt(end).*1e-3);
        
        [trlID] = sort([lfpDat.hitIdx;lfpDat.missIdx]);
        [hitIdx] = find(ismember(lfpDat.hitIdx,trlID));
        [missIdx] = find(ismember(lfpDat.missIdx,trlID));
        [ntrl] = length(lfpDat.trlENC);
        
        if ~(length(trlID) == size(lfp,2))
            error('trial selection is out of range');
        end;
        
        if length(dt) ~= size(lfp,1)
            error('matrix dimensions out of range');
        end;
        
        if ~(length(trlID) == ntrl) ||  ~( size(lfp,2) == ntrl )
            error('trial selection is out of range');
        end;
        
        %%
        bsz=(dt2(4)-dt2(3))*1e-3;
        
        ixB = [find(dt>=-900 & dt <=-650);...
            find(dt>=-600 & dt <=-350);...
            find(dt>=-300 & dt <=-50)];
        
        lag = 500;
        
        nTw1 = length(250:25:4000);
        nTw2 = length(2000:25:4000);
        
        cnt = 0;
        [spkSelIdx] = [];
        [ mxTE ] = [];
        [ mxDE ] = [];
        [spkRaster] = {};
        [xc] = {};
        [wvf] = {};
        [isiH] = {};
        [n2] = {};
        
        for kt = 1:length(lfpDat.chanLab)
            
            spk = spkDat.sortedSpikes{kt}
            
            spkCnt = spk.SpikeTimesSeg.*1e3;
            spkCnt = spkCnt(spkCnt>=0 & spkCnt<=4000);
            
            if length(spkCnt) >=50
                
                cnt = cnt+1;
                spkSelIdx(cnt) = kt;
                tsOrig = (spk.newSpikeTimes./1e6).*1e3;
                dx = diff(tsOrig);
                
                [isiH{cnt}] = hist(dx,0:501);
                isiH{cnt}(end) = [];
                
                wvf{cnt} = spk.wavf(spk.oriIx,:);
                
                spkRaster{cnt} = zeros(length(trlID),length(dt));
                n2{cnt} = zeros(length(trlID),length(dt2));
                for it = 1:length(trlID)
                    ts = spk.SpikeTimesSeg(spk.trl == trlID(it));
                    ts = ts.*1e3;
                    n2{cnt}(it,:) = hist(ts,dt2);
                    ts = ts(ts>=dt(1) & ts <=dt(end));
                    spkRaster{cnt}(it,:) = hist(ts,dt);
                end;
                n2{cnt}(:,[1 end]) = [];
                
                xc{cnt} = zeros(size(spkRaster{cnt},1),lag*2+1);
                for jt = 1:size(spkRaster{cnt},1)
                    xc{cnt}(jt,:) = xcorr(spkRaster{cnt}(jt,:),lag);
                end;
                xc{cnt}(:,lag+2:end);
                
                [baseSpkCnt] = zeros(length(trlID),size(ixB,1));
                for it = 1:size(ixB,1)
                    [baseSpkCnt(:,it)] = sum(spkRaster{cnt}(:,ixB(it,:)),2);
                end;
                [baseSpkCnt] = median(baseSpkCnt,2);
                
                [ix1] = find(dt>=0 & dt <=250);
                
                empT = zeros(1,nTw1);
                empD = zeros(1,nTw2);
                tt2 = 0;
                for tt = 1:nTw1
                    if tt>1
                        ix1 = ix1+25;
                    end;
                    [postSpkCnt] = sum(spkRaster{cnt}(:,ix1),2);
                    if sum(postSpkCnt)>=fix(length(trlID)/4)
                        [~,~,~,stats] = ttest(postSpkCnt,(baseSpkCnt));
                        empT(tt) = stats.tstat;
                    else
                        empT(tt) = 0;
                    end;
                    if dt(ix1(1))>=2000
                        tt2 = tt2+1;
                        for st = 1:length(stimID)
                            ix = find(strcmp(stimCat , stimID(st)));
                            empD(tt2) = empD(tt2)+sum(postSpkCnt(ix).*(log(mean(postSpkCnt(ix)))-log(mean(postSpkCnt))));
                        end;
                        empD(tt2) =  2.*empD(tt2);
                    end;
                    
                end;
                mxTE(cnt) = max(empT);%empT(min(find(abs(empT) == max(abs(empT)))));
                mxDE(cnt) = max(empD);
            end;
            
        end;
        
        %%
        if length(spkSelIdx) ~= length(spkRaster)
            error('data must have matching dimensions');
        end;
        
        [randT] = zeros(nrand,length(spkSelIdx));
        [randD] = zeros(nrand,length(spkSelIdx));
        tic;
        parfor it = 1:nrand
            fprintf([num2str(it),'/',num2str(nrand)]);
            rT = zeros(length(spkSelIdx),1);
            rD = zeros(length(spkSelIdx),1);
            rStimCat = stimCat(randperm(length(stimCat)));
            for kt = 1:length(spkSelIdx)
                
                [ix1] = find(dt>=0 & dt <250);
                
                tstat = zeros(1,nTw1);
                rd  = zeros(1,nTw2);
                tt2 = 0;
                for tt = 1:nTw1%:7%:6
                    if tt>1
                        ix1 = ix1+25;
                    end;
                    
                    x1 =[sum(spkRaster{kt}(:,ixB),2) sum(spkRaster{kt}(:,ix1),2)];
                    
                    ixa = 1:ntrl;
                    ixb = ntrl+1:ntrl*2;
                    
                    x1 = x1(:);
                    temp1 = x1(ixa);
                    temp2 = x1(ixb);
                    rIx = randperm(ntrl);
                    rIx = rIx(1:2:end);
                    x1(ixa(rIx)) = temp2(rIx);
                    x1(ixb(rIx)) = temp1(rIx);
                    x1 = reshape(x1,[ntrl 2]);
                    
                    [~,~,~,stats] = ttest(x1(:,2),x1(:,1));
                    tstat(tt) =  stats.tstat;
                    
                    if dt(ix1(1))>=2000
                        tt2 = tt2+1;
                        x1 = sum(spkRaster{kt}(:,ix1),2);
                        for st = 1:length(stimID)
                            ix = find(strcmp(rStimCat,stimID(st)));
                            rd(tt2) = rd(tt2)+ sum(x1(ix).*(log(mean(x1(ix)))-log(mean(x1(:)))));
                        end;
                        rd(tt2) = 2.*rd(tt2);
                    end;
                end;
                
                rT(kt) = max(tstat);
                rD(kt) = max(rd);
            end;
            randT(it,:) = rT;
            randD(it,:) = rD;
            fprintf('\n');
        end;
        
        %%
        pval1 = zeros(length(mxTE),1);
        for it = 1:length(mxTE)
            pval1(it) = length(find(randT(:,it)>=mxTE(it)))/size(randT,1);
        end;
        
        pval2 = zeros(length(mxDE),1);
        for it = 1:length(mxDE)
            pval2(it) = length(find(randD(:,it)>=mxDE(it)))/size(randD,1);
        end;
        toc;
        
        %%
        % for kt = 1:length(spkSelIdx)
        %
        %     if  (pval1(kt)<=5e-2) || (pval2(kt)<=5e-2) %|| (pval2(kt)<=5e-2)
        %
        %         n = spkRaster{kt};
        %
        %         figure;
        %         subplot(221);
        %         plot(0:500,isiH{kt});xlim([-10 500]);
        %         title(lfpDat.chanLab{spkSelIdx(kt)});
        %         subplot(222);
        %         hold on;
        %         for jt = 1:size(n,1)
        %             y = jt*ones(1,sum(n(jt,:)==1));
        %             x = dt(n(jt,:)==1);
        %             x = [x;x];
        %             y = [y-.5;y+.5];
        %             line(x,y,'Color','k');
        %         end;
        %         axis tight;
        %         plot([0 0],[0 size(n,1)+1],'r');
        %         plot([2 2].*1e3,[0 size(n,1)+1],'r');
        %         xlim([dt(1) dt(end)]);
        %         subplot(224);
        %         hold on;
        %         fr = sum(n2{kt},1)./size(n2{kt},1)./bsz;
        %         plot(dt2(2:end-1),fr,'s-','Color',[.75 .75 .75],'MarkerFaceColor','c');
        %         plot([0 0],[min(fr) max(fr)],'k');
        %         plot([2 2].*1e3,[min(fr) max(fr)],'k');
        %         axis tight;xlim([dt2(2) dt2(end-1)]);
        %         subplot(223);
        %         plot(1:lag,mean(xc{kt},1));
        %         axis tight;
        %     end;
        % end;
        
        %%
        sigSel = find(sum([pval1 pval2]<0.05,2)~=0);
        fr = zeros(length(sigSel),length(dt2)-2);
        for it = 1:length(sigSel)
            
            x = sum(n2{sigSel(it)},1)./size(n2{sigSel(it)},1)./bsz;
            x = (x-mean(x(dt2<=0)))./std(x(dt2<=0));
            fr(it,:) = x;
        end;
        
        %%
        savePath = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',pID{pt},'/spkRes/'];
        saveName = [pID{pt},'_',expMode,'_',seshID{seshIt},'_spkResp.mat'];
        
        [hitIdx] = lfpDat.hitIdx;
        [missIdx] = lfpDat.missIdx;
        
        save([savePath,saveName],'sigSel','spkSelIdx','spkRaster','isiH','n2','pval1','pval2','mxTE','mxDE','randT','randD','fr','wvf','missIdx','hitIdx');
        
        % %%
        % M = mean(fr,1);
        % SE = std(fr,0,1)./sqrt(size(fr,1)-1);
        %
        % figure;
        % jbfill(dt2(2:end-1),M-SE,M+SE,[.9 0 0],[.9 0 0],1,.4);
        % hold on;
        % plot(dt2(2:end-1),M,'r','LineWidth',3);
        % xlim([-250 dt2(end)]);
    end;
end;

%%
delete(gcp);
clear;
exit;
