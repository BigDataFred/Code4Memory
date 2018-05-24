function [MWdat] = localRefMWdat(MWdat)
rms = [];
for jt = 1:length( MWdat.trial )
    
    rms(jt,:) = sqrt(sum(MWdat.trial{jt}.^2,2));
    
end;

rms = mean(rms,1);

figure;
idx = 1:8;
ixRef = [];
c = 0;
for it = 1:8:size(rms,2)
    c = c+1;
    
    [~,mIx] = min(rms);
    
    ixRef(c) = idx(mIx);
    
    hold on;
    plot(idx,rms,'s-','Color','k');
    plot(idx(mIx),rms(mIx),'ys','MarkerFaceColor','y');
    
    idx = idx+8;
    
end;
x = MWdat.label;
for it = 1:length(x);x{it}(regexp(x{it},' ')) = [];end;
for it = 1:length(x);x{it}(regexp(x{it},'_')) = [];end;

ylabel('\Delta-RMS');
set(gca,'XTick',2:2:length(MWdat.label));
set(gca,'XTickLabel',MWdat.label(2:2:end));
set(gca,'XTickLabelRotation',-45);

cfg                     = [];
cfg.channel             = MWdat.label(ixRef);

[refDat] = ft_selectdata( cfg , MWdat );

cfg                     = [];
cfg.channel             = setdiff(1:length(MWdat.label),ixRef);

[MWdat] = ft_selectdata( cfg , MWdat );


for it = 1:length(MWdat.trial)
    MWdat.trial{it} = MWdat.trial{it} - repmat(refDat.trial{it},[length(MWdat.label) 1]);
end;