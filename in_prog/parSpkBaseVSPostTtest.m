%%
if isempty(gcp('nocreate'))
    parpool(36,'SpmdEnabled',false);
end;

[trlID] = sort([lfpDat.hitIdx;lfpDat.missIdx]);

dt = [-1e3:5e3];

[tw1] = createVariableTimeWindowsAcrossBinRange([-1000 0],10,20);
[tw2] = createVariableTimeWindowsAcrossBinRange([3000 4000],10,20);

mxTE = [];
spkRaster = {};
cnt = 0;
selIx = [];
for kt = 1:length(lfpDat.chanLab)
        
    spkTs = spkDat.sortedSpikes{kt}.SpikeTimesSeg.*1e3;    
    if sum(spkTs>=dt(1) & spkTs<=dt(end)) >=100
        
        cnt = cnt+1;                
        selIx(cnt) = kt;
        spkRaster{cnt} = zeros(length(trlID),length(dt));
        for it = 1:length(trlID)
            ts = spkDat.sortedSpikes{kt}.SpikeTimesSeg(spkDat.sortedSpikes{kt}.trl == trlID(it));
            ts = ts.*1e3;
            ts = ts(ts>=dt(1) & ts <=dt(end));
            spkRaster{cnt}(it,:) = hist(ts,dt);
        end;               
        
        tEMP = zeros(size(tw1,1),1);
        parfor it = 1:size(tw1,1)
            ixB = find(dt>= tw1(it,1) & dt <= tw1(it,2));
            ixP = find(dt>= tw2(it,1) & dt <= tw2(it,2));
            spkCountB = sum(spkRaster{cnt}(:,ixB),2);
            spkCountP = sum(spkRaster{cnt}(:,ixP),2);
            [~,~,~,stat] = ttest(spkCountP,spkCountB);
            tEMP(it) = stat.tstat;
        end;
        
        mxTE(cnt,2) = [min(tEMP) max(tEMP)];% maximum across time windows per channel
        
    end;
end;

%%
[ntrl] = size(spkRaster{1},1);

[twIx] = zeros(size(tw,1),length(selIx));
for it = 1:size(tw1,1)
    x =  find(dt>= tw1(it,1) & dt <= tw1(it,2));
    twIx(it,1:2) = [min(x) max(x)];
    x =  find(dt>= tw2(it,1) & dt <= tw2(it,2));
    twIx(it,3:4) = [min(x) max(x)];
end;

cnt = 0;
x = zeros(ntrl*length(selIx)*size(tw1,1),2);
ix = 1:ntrl;
for kt = 1:length(selIx)
    n = spkRaster{kt};
    
    for it = 1:size(tw1,1)
        x(ix,:) = [sum(n(:,twIx(it,1):twIx(it,2)),2) sum(n(:,twIx(it,3):twIx(it,4)),2)];
        ix = ix+ntrl;
    end;
end;

tic;
tEMPrd = zeros(1e4,2);
parfor jt = 1:1e4
    tRD = zeros(size(x,1)/ntrl,1);
    ix = 1:ntrl;
    for it = 1:size(x,1)/ntrl
        spk = x(ix,:);
        spkr = spk(:);
        ix1 = 1:ntrl;% column1
        ix2 = ntrl+1:ntrl*2;% column2  
        
%         if any(spkr(ix1)~=spk(:,1)) || any(spkr(ix2) ~= spk(:,2))
%             error('indexes are out of range');
%         end;
        
        [temp1] = spkr(ix1);%base
        [temp2] = spkr(ix2);%post
        
        rIx = randperm(ntrl);% shuffle indexes
        rIx = rIx(1:2:end);% 50% of trials
        
        spkr(ix1(rIx)) = temp2(rIx);% randomly flip bl and post
        spkr(ix2(rIx)) = temp1(rIx);
        
        spkr = reshape(spkr,[ntrl 2]);                      
        spkCountB = spkr(:,1);
        spkCountP = spkr(:,2);
        [~,~,~,stat] = ttest(spkCountP,spkCountB);
        tRD(it) = stat.tstat;
        ix = ix+ntrl;
    end;
    tEMPrd(jt,2) = [min(tRD) max(tRD)];% maximum across time windows and channels
end;
toc;

%%
pval1 = zeros(length(selIx),1);
pval2 = zeros(length(selIx),1);
for it = 1:length(selIx)
    pval1(it) = length(find(tEMPrd(:,2)>=mxTE(it,2)))/1e4;
    %pval2(it) = length(find(tEMPrd(:,1)<=mxTE(it,1)))/1e4;
end;



