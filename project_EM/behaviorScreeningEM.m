%%
pID = {'P02','P04','P05','P07','P08'};

bL  = cell(length(pID),1);
pCR = cell(length(pID),1);
tpR  = cell(length(pID),1);
fpR = cell(length(pID),1);
H1 = [];RT1 = [];
H2 = [];RT2 = [];
H3 = [];RT3 = [];
for currPat = 1:length(pID)
    h1 = []; h2 = []; h3 = [];
    rt1 = []; rt2 = []; rt3 = [];
    expMode = 'fVSpEM';
    
    chck = dir(['/media/rouxf/rds-share/iEEG_DATA/MICRO/',pID{currPat},'/',expMode,'/']);
    chck(1:2) = [];
    cnt = 0; sel = [];
    for it = 1:length(chck)
        if  ( regexp(chck(it).name,'\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}') ) && ( chck(it).isdir )
            cnt = cnt+1;
            sel(cnt) = it;
        end;
    end;
    sesh = {chck(sel).name}'
       
    for seshIt = 1:length(sesh)
        p2d = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',pID{currPat},'/',expMode,'/',sesh{seshIt},'/'];
        fn = dir([p2d,pID{currPat},'_',sesh{seshIt},'_LogFile_EMtask_LogDat.mat']);
        
        load([p2d,fn.name])
        
        stim = logDat.LogDat1.log(:,3:4);
        resps = str2double(logDat.LogDat2.log(:,5:6));
        rts = str2double(logDat.LogDat2.log(:,8:9));

        lab = {};
        for jt = 1:size(stim,1)
            lab(jt,:) = {stim{jt,1}(1) stim{jt,2}(1)};
        end;
        
        ix1 = [];ix2 = [];ix3 = [];
        for jt = 1:size(lab,1)
            if strcmp(lab(jt,:),{'f' 'f'});
                ix1 = [ix1 jt];
            elseif strcmp(lab(jt,:),{'p' 'p'});
                ix2 = [ix2 jt];
            else
                ix3 = [ix3 jt];
            end;
        end;
        x =resps(ix1,:);
        x = x(:);
        h1 = [h1 sum(x)./length(x)];
        rt1 = [rt1 mean(median(rts(ix1,:)))];
        
        x =resps(ix2,:);
        x = x(:);
        h2 = [h2 sum(x)./length(x)];
        rt2 = [rt2 mean(median(rts(ix2,:)))];
        
        x =resps(ix3,:);
        x = x(:);
        h3 = [h3 sum(x)./length(x)];
        rt3 = [rt3 mean(median(rts(ix3,:)))];
        
        bL{currPat} = [bL{currPat};diff(logDat.LogDat1.idx,[],2)+1];
        
        x = str2double(logDat.LogDat2.log(:,5:6));
        for it = 1:size(logDat.LogDat2.idx,1)
            ix = logDat.LogDat2.idx(it,1):logDat.LogDat2.idx(it,2);
            sel = x(ix,:);
            RTs = str2double(logDat.LogDat2.log(ix,8:9));
            delIx = find(mean(RTs,2) ==1e6);
            sel(delIx,:) = [];
            tpR{currPat} = [tpR{currPat} length(find(sel==1))/length(sel(:))];
            fpR{currPat} = [fpR{currPat} length(find(sel==0))/length(sel(:))];
            
            p1 = (length(find(sel(:,1)==1))/size(sel,1)-(1/4))/(1-(1/4));
            p2 = (length(find(sel(:,2)==1))/size(sel,1)-(1/3))/(1-(1/3));
            pCR{currPat} = [pCR{currPat} p1*p2];
        end;
        
    end;
    H1 = [H1 mean(h1)];
    H2 = [H2 mean(h2)];
    H3 = [H3 mean(h3)];
    RT1 = [RT1 mean(rt1)];
    RT2 = [RT2 mean(rt2)];
    RT3 = [RT3 mean(rt3)];
end;

%%
cnt = 0;
figure;
for curPat = 1:length(pID)
    cnt = cnt+1;    
    subplot(length(pID),3,cnt);
    hold on;
    plot(1:length(bL{curPat}),bL{curPat},'LineWidth',3);
    xlabel('Learning block');
    ylabel('# of pairs');
    xlim([0 length(bL{curPat})+1]);
    
    cnt = cnt+1;
    subplot(length(pID),3,cnt);
    hold on;
    plot(1:length(pCR{curPat}),pCR{curPat},'LineWidth',3);
    plot([1 length(pCR{curPat})],[1/4*1/3 1/4*1/3],'r--');
    xlabel('Learning block');
    ylabel('PCR');
    xlim([0 length(bL{curPat})+1]);
    ylim([0 1]);
    
    cnt = cnt+1;
    subplot(length(pID),3,cnt);
    x = jet(length(bL{curPat}));
    hold on;
    [~,sIx] = sort(bL{curPat});
    for it = 1:length(sIx)
        plot(fpR{curPat}(sIx(it)),tpR{curPat}(sIx(it)),'s','Color',x(it,:));
    end;
    plot([0 1],[0 1],'k');
    ylabel('TPR');
    xlabel('FPR');
end;

figure;
subplot(121);
hold on;
x = [H1' H2' H3'];
for it = 1:size(x,2)
    plot(it*ones(size(x(:,it),1),1),x(:,it),'ow','MarkerFaceColor','b','MarkerSize',8);
end;
xlim([0 4]);
ylim([0 1]);
set(gca,'XTick',[1:3]);
set(gca,'XTickLabel',{'ff' 'pp' 'fp'});
ylabel('Percent CR (%)');
xlabel('Stimulus pair');
subplot(122);
hold on;
x = [RT1' RT2' RT3'];
for it = 1:size(x,2)
    plot(it*ones(size(x(:,it),1),1),x(:,it),'ow','MarkerFaceColor','b','MarkerSize',8);
end;
xlim([0 4]);
ylim([0 5]);
set(gca,'XTick',[1:3]);
set(gca,'XTickLabel',{'ff' 'pp' 'fp'});
ylabel('Reaction time (s)');
xlabel('Stimulus pair');