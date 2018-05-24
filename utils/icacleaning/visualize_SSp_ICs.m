function visualize_SSp_ICs(SSpidx,protoBlink,bcComp,comp_data)
if ~isempty(SSpidx)
    %%
    if mod(length(SSpidx),2) == 0
        n = length(SSpidx);
    else
        n = length(SSpidx)+1;
    end;
    figure;
    set(gcf,'Position',[900 5 950 1050]);
    idx1 = zeros(1,ceil(length(SSpidx)/2));for it = 1:ceil(length(SSpidx)/2);if it <2;idx1(it) = it;else idx1(it) = idx1(it-1)+4;end;end;
    subplot(n,4,[idx1]);
    a = gca;
    imagesc(linspace(-.5,.5,201),1:size(protoBlink{1},1),squeeze(protoBlink{1}));
    
    idx3 = idx1(end)+4:4:idx1(end)+4+(4*(ceil(length(SSpidx)/2)-1));
    subplot(n,4,[idx3]);
    a = [a gca];
    imagesc(linspace(-.5,.5,201),1:size(protoBlink{2},1),squeeze(protoBlink{2}));
    
    
    idx2 = zeros(1,ceil(length(SSpidx)/2));for it = 1:ceil(length(SSpidx)/2);if it <2;idx2(it) = it+1;else idx2(it) = idx2(it-1)+4;end;end;
    subplot(n,4,[idx2]);
    a = [a gca];
    hold on;
    plot([0 0],[min(squeeze(mean(protoBlink{1},1))) max(squeeze(mean(protoBlink{1},1)))],'r--');
    plot(linspace(-.5,.5,201),squeeze(mean(protoBlink{1},1)),'k','LineWidth',3);axis tight;
    title(['Peak triggered average vEOG']);
    
    idx4 = idx2(end)+4:4:idx2(end)+4+(4*(ceil(length(SSpidx)/2)-1));
    subplot(n,4,[idx4]);
    a = [a gca];
    hold on;
    plot([0 0],[min(squeeze(mean(protoBlink{2},1))) max(squeeze(mean(protoBlink{2},1)))],'r--');
    plot(linspace(-.5,.5,201),squeeze(mean(protoBlink{2},1)),'k','LineWidth',3);axis tight;
    title(['Peak triggered average hEOG']);
    
    for it = 1:length(a)
        xlabel(a(it),'Time [s]');
    end;
    for it = 1:length(a)-2
        ylabel(a(it),'Trial #');
    end;
    for it = 3:length(a)
        ylabel(a(it),'Amplitude [\sigma]');
    end;
    %%
    cfg = [];
    cfg.layout = 'neuromag306planar.lay';
    cfg.comment = 'no';
    cfg.marker = 'off';
    
    k = setdiff(1:(ceil(length(SSpidx))*4),sort([idx1 idx2 idx3 idx4]));
    
    subplot(n,4,k(1:2:end));
    hold on;
    a = gca;
    hold on;
    plot([0 0],[min(min(squeeze(mean(bcComp{2}(SSpidx,:,:),2)))) max(max(squeeze(mean(bcComp{2}(SSpidx,:,:),2))))],'r--');
    plot(linspace(-.5,.5,201),squeeze(mean(bcComp{2}(SSpidx,:,:),2)),'k');
    axis tight;
    
    c=0;
    for it = 2:2:length(k);
        subplot(n,4,k(it));
        c = c+1;
        cfg.component = SSpidx(c);
        ft_topoplotIC(cfg,comp_data);
    end;
    
    for it = 1:length(a);
        xlabel(a(it),'Time [s]');
        ylabel(a(it),'Amplitude [\sigma]');
    end;
end;

