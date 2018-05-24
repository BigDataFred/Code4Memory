function [D1,D2] = stimulusSelectivityIndex(spkDat,stimIdx)
%%
b = 100:10:990;
%%
c =0;
for it = 1:length(b)
    for jt = it+1:length(b)
        c=c+1;
        p(c,:) = [it jt];
    end;
end;

%%
nperms = 1e3;

%%
D1 = cell(1,size(p,1));
D2 = cell(size(p,1),nperms);

nIx = []; for jt = 1:length(stimIdx);nIx(jt) = length(stimIdx{jt});end;

parfor it = 1:size(p,1)
    
    dt1 = b(p(it,1))-10:b(p(it,1))+10;
    dt2 = b(p(it,2))-10:b(p(it,2))+10;
    
    scTot = [];
    for jt = 1:length(spkDat)
        spkCount = length(find(spkDat{jt}>=dt2(1) & spkDat{jt}<=dt2(end)));
        scTot(jt) = spkCount;%total spike count 
    end;
    
    for n = 1:length(stimIdx)
        
        scStim = [];
        for jt = 1:length(stimIdx{n})
            spkCount = length(find(spkDat{stimIdx{n}(jt)} >= dt1(1) & spkDat{stimIdx{n}(jt)} <= dt1(end)));
            scStim(jt) = spkCount;% sipke count for specific stimulus     
        end;
        
        D1{it}(n) = nansum(scStim.*(log(mean(scStim))-log(mean(scTot))));
        
    end;
    D1{it} = 2*nansum(D1{it});

    
    dum = cell(1,nperms);
    for jt = 1:nperms
        
        stimIdx2 = [];
        ixPool = randperm(length(spkDat));
        stimIdx2 = {ixPool(1:nIx(1)),ixPool(nIx(1)+1:nIx(1)+nIx(2)),ixPool(nIx(1)+nIx(2)+1:nIx(1)+nIx(2)+nIx(3))};
        
        for n = 1:length(stimIdx2)
            
            scStim = [];
            for kt = 1:length(stimIdx2{n})
                spkCount = length(find(spkDat{stimIdx2{n}(kt)} >= dt1(1) & spkDat{stimIdx2{n}(kt)} <= dt1(end)));
                scStim(kt) = spkCount;% sipke count for specific stimulus
            end;
            
            dum{jt}(n) = nansum(scStim.*(log(mean(scStim))-log(mean(scTot))));
            
        end;
        dum{jt} = 2*nansum( dum{jt});
    end;
    D2(it,:) = dum(:);
end;

return;