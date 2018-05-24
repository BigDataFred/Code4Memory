npow = {};
for it = 1:length(pow)
    
    cfg                 =  [];
    cfg.baseline         = [-0.5 0];
    cfg.baselinetype    = 'absolute';
    
    npow{it}             = ft_freqbaseline(cfg , pow{it} );
    
end;

%%
m1 = min(min(min(min(npow{1}.powspctrm))));
m2 = max(max(max(min(npow{1}.powspctrm))));

m1 = m1/100*1;
m2 = m2/100*1;

n = size(pow{1}.powspctrm,2)/8;

figure;
for it = 1:length(npow{1}.label)
    subplot(n,8,it);
    hold on;
    pcolor(npow{1}.time,npow{1}.freq,squeeze(mean(log10(pow{1}.powspctrm(ix{4},it,:,:)),1)-mean(log10(pow{1}.powspctrm(ix{6},it,:,:)),1)));shading interp;
    plot([-1 -1],[min(npow{1}.freq) max(npow{1}.freq)],'w');
    plot([0 0],[min(npow{1}.freq) max(npow{1}.freq)],'w');
    plot([2 2],[min(npow{1}.freq) max(npow{1}.freq)],'w');
    axis tight;
    %caxis([m1 m2]);
end;
%%

n = size(pow{2}.powspctrm,2)/8;

figure;
for it = 1:length(npow{2}.label)
    subplot(n,8,it);
    hold on;
    pcolor(npow{2}.time,npow{2}.freq,squeeze(mean(pow{2}.powspctrm(ix{4},it,:,:),1)-mean(pow{2}.powspctrm(ix{6},it,:,:),1)));shading interp;
    plot([-1 -1],[min(npow{2}.freq) max(npow{2}.freq)],'w');
    plot([0 0],[min(npow{2}.freq) max(npow{2}.freq)],'w');
    plot([2 2],[min(npow{2}.freq) max(npow{2}.freq)],'w');
    axis tight;
    %caxis([m1 m2]);
end;