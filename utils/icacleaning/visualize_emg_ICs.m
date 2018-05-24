function visualize_emg_ICs(comp_freq,comp,emg_idx)
%%
figure;
set(gcf,'Position',[900 5 950 1050]);
k = 0;
if ~isempty(emg_idx)
    
    for it = 1:length(emg_idx)
        k = k+1;
        
        subplot(length(emg_idx),2,k); 
        hold on;
        plot(comp_freq.freq,squeeze(mean(log10(comp_freq.powspctrm),1)),'Color',[.75 .75 .75]);
        plot(comp_freq.freq,squeeze(mean(log10(comp_freq.powspctrm(:,[emg_idx(k)],:)),1)),'r');
        axis tight;
        xlabel('Frequency [Hz]');
        ylabel('Power [log10]');
        
        subplot(length(emg_idx),2,k+1);                
        cfg = [];
        cfg.layout = 'neuromag306planar.lay';        
        cfg.comment = 'no';
        cfg.electrodes = 'off';
        cfg.component = emg_idx(k);
        ft_topoplotIC(cfg,comp);
        title(['component #',num2str(emg_idx(k))]);
    end;
    
%     %%
%     dum = struct;
%     dum.cfg = [];
%     dum.label = comp.label(emg_idx);
%     dum.trial = cell(1,length(comp.trial));
%     dum.time = cell(1,length(comp.time));
%     for jt = 1:length(comp.trial)
%         dum.trial{jt} = (comp.trial{jt}(emg_idx,:));
%         dum.time{jt} = comp.time{jt};
%     end;
%     cfg = [];
%     cfg.derivative = 'yes';
%     
%     dum = ft_preprocessing(cfg,dum);
%     
%     cfg = [];
%     cfg.method = 'mtmconvol';
%     cfg.foi = 20:100;
%     cfg.taper = 'dpss';
%     cfg.padding = 'maxperlen';
%     cfg.tapsmofrq = 15;
%     cfg.t_ftimwin = 0.25*ones(1,length(cfg.foi));
%     cfg.toi = dum.time{1}(1):0.025:dum.time{1}(end);
%     
%     emg_freq = ft_freqanalysis(cfg,dum);   
%     %%
%     for jt = 1:size(emg_freq.powspctrm,1)
%         subplot(2,length(emg_idx)+3,length(emg_idx)+3+3+jt);
%         imagesc(emg_freq.time,emg_freq.freq,(squeeze((emg_freq.powspctrm(jt,:,:)))));
%         axis xy;
%         xlabel('Time [s]');
%         ylabel('Frequency [Hz]');
%         title(['Component #',emg_freq.label{jt}(end-2:end)]);
%         
%     end;
end;
