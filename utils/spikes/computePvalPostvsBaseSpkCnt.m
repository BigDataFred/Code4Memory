function [pval] = computePvalPostvsBaseSpkCnt(spkResp) 
%%
dt = -1500:5500;
dt2 = -1500:250:5500;

ntrl = length( spkResp.spikeTimes );

sigSel = [];

[spkRaster] = zeros(ntrl,length(dt));
[iFR] = zeros(ntrl,length(dt));
[fr] = zeros(ntrl,length(dt2));
for curTrl = 1:ntrl
    x = spkResp.spikeTimes{curTrl};
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

iFR = mean(iFR,1);

xix1 = find(sign(diff(iFR(tx>=200)>trsh))==1)+1;
xix2 = find(sign(diff(iFR(tx>=200)>trsh))==-1)+1;

if ~isempty(xix1) || ~isempty(xix2)
    if min(xix2) < min(xix1)
        xix1 = [1 xix1];
    end;
    if length(xix1)>length(xix2)
        xix2(end+1) = length(iFR);
    end;
end;

[baseCnt] = sum(spkRaster(:,tx>=-400 & tx <=0),2);
baseSD = 2*std(baseCnt);

baseLine = median(baseCnt);

pval = 1;
if any(diff([xix1' xix2'],[],2)>100)
    
    TW = [200 600];
    p = zeros(1,44);
    h = zeros(1,44);
    for curTW = 1:44
        [postCnt] = sum(spkRaster(:,tx>=TW(1) & tx <=TW(2)),2);
        
        [p(curTW)] = ranksum(postCnt(:),baseCnt(:));
        TW = TW+100;
        h(curTW) = 1;
    end;
    pval = min(p(h==1));
    
   
end;

