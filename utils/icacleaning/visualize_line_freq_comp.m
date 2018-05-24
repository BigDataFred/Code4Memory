function visualize_line_freq_comp(comp,comp_freq,line_idx)

%%
if line_idx >0
    cfg = [];
    cfg.layout = 'neuromag306planar.lay';
    cfg.component = [];
    cfg.marker = 'off';
    cfg.comment = 'no';
    
    figure;
    set(gcf,'Position',[900 5 950 1050]);
    
    idx = zeros(1,ceil(length(line_idx)/4));
    for jt = 1:ceil(length(line_idx)/4)
        idx(jt) = (jt-1)*ceil(length(line_idx)/4)+1;
    end;
    
    idx = [idx idx+1];
    if length(line_idx > 1)
        subplot(ceil(length(line_idx)/4),6,idx);
    else
        subplot(1,6,idx);
    end;
    hold on;
    plot(comp_freq.freq,squeeze(mean((comp_freq.powspctrm),1)),'Color',[.75 .75 .75]);
    h = [];
    h = plot(comp_freq.freq,squeeze(mean(mean((comp_freq.powspctrm(:,line_idx,:)),2),1)),'r');
    xlabel('Frequency [Hz]');
    ylabel('Power [log10]');
    axis tight;
    legend(h,num2str(line_idx));
    
    idx = setdiff(1:ceil(length(line_idx)/4)*6,idx);
    for it = 1:length(line_idx)
        cfg.component = line_idx(it);
        subplot(ceil(length(line_idx)/4),6,idx(it));
        ft_topoplotIC(cfg,comp);
        title(['IC #',num2str(line_idx(it))]);
    end;
    
else;
    figure;
    clf;
end;