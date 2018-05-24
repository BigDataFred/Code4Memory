function [LFPsig,chID] = extractMedianLFP(dat,LFPr,chID)

if isstruct(dat)
    f =1;
    LFPsig = dat;
    LFPsig.label = {};
    LFPsig.trial = {};
    ntrl = length(dat.trial);
    nsamp = length(dat.trial{1});
else
    f = 2;
    [nsamp,ntrl] = size(dat{1});
    LFPsig = cell(1,1);
    x = zeros(nsamp,length(chID),ntrl);
    for it = 1:length(chID)
        x(:,it,:) = dat{it};
    end;
    
end;

chLab = chID;
chID = unique(chID);

selIx = cell(1,length(chID));
for it = 1:length(chID)
    chIdx = find(strcmp(chLab,chID(it)));
    dum = squeeze(mean(LFPr(:,chIdx,chIdx),1));
    dum = dum - eye(size(dum));
    dum(sign(dum)==-1) = 0;
    dum = dum >.7;
    selIx{it} = find(sum(dum,2) == mode(sum(dum,2)));
    if f==1
        LFPsig.label(it) = {['median',chID{it}]};
    end;
end;

for jt = 1:length(chID)
        
    for it = 1:ntrl
        if f==1
            LFPsig.trial{it}(jt,:) = median(dat.trial{it}(selIx{jt},:),1);
        elseif f==2
            LFPsig{jt}(:,it) = median(x(:,selIx{jt},it),2);
        end;
        if any(isnan(x))
            error('timesries contains NaNs');
        end;
    end;
    
end;
