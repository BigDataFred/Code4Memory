%%
cfg = [];
cfg.channel = 'all';%comp.label(comp_idx4);
cfg.keeptrials = 'yes';

SD_EP = ft_timelockanalysis(cfg,comp);
%%
figure;
subplot(2,4,1:3);
imagesc(SD_EP.time,1:size(SD_EP.trial,1),(squeeze(SD_EP.trial(:,s_idx(min(find(sd==min(sd)))),:))));
xlabel('Time [s]');
ylabel('Trial #');
axis tight;
axis xy;
title(['Example of stationary IC']);

subplot(2,4,5:7);
imagesc(SD_EP.time,1:size(SD_EP.trial,1),(squeeze(SD_EP.trial(:,s_idx(max(find(sd==max(sd)))),:))));
xlabel('Time [s]');
ylabel('Trial #');
axis xy;
axis tight;
title(['Example of non-stationary IC']);

subplot(244);
plot((SD(s_idx(find(sd==min(sd))),:)),1:size(SD,2));
axis tight;
set(gca,'XTick',[]);
%xlim([min(SD(find(sd==max(sd)),:)) max(SD(s_idx(find(sd==max(sd))),:))]);
title('Variance across samples');

subplot(248);
plot((SD(s_idx(find(sd==max(sd))),:)),1:size(SD,2));
axis tight;
set(gca,'XTick',[]);
%xlim([min(SD(find(sd==max(sd)),:)) max(SD(s_idx(find(sd==max(sd))),:))]);
title('Variance across samples');
%     saveas(gcf,['/data/home1/froux/DATASETS/DevelopmentData/figures/ICA_cleaning/',files1(ot).name(1:5),'_comp_idx4_1.fig'],'fig');
%     close;
%     %%
%     figure;
%     set(gcf,'Position',[50 130 1730 820]);
%     subplot(221);
%     hold on;
%     plot(1:length(sd),sd,'ks','MarkerFaceColor','k');axis tight;
%     plot([min(idx) min(idx)],[min(sd) max(sd)],'r--');
%     plot([1 length(sd)],[2 2],'r--');
%     ylabel('Variance [z-score]');
%     xlabel('IC # sorted by Variance');
%
%     subplot(223);
%     hold on;
%     plot(1:length(roc),roc,'-','Color','g');axis tight;
%     plot([min(idx) min(idx)],[min(roc) max(roc)],'r--');
%     plot([1 length(roc)],[2 2],'r--');
%     ylabel('Rate of change [z-score]');
%     xlabel('IC # sorted by Variance');
%
%     subplot(222);
%     hold on;
%     dum = SD;
%     dum = dum(s_idx,:);
%     imagesc(log10(dum));
%     axis xy;
%     if ~isempty(comp_idx4)
%         plot([1 size(SD,2)],[size(SD,1)-length(comp_idx4) size(SD,1)-length(comp_idx4)],'w--','LineWidth',3);
%     end;
%     axis tight;
%     xlabel('Trial #');
%     ylabel('IC # sorted by Variance');
%     title('Variance of IC over time');
%     colormap jet;
%
%     subplot(224);
%     hold on;
%     plot(comp_freq.freq,squeeze(mean(20.*log10(comp_freq.powspctrm),1)),'Color',[.75 .75 .75]);
%     if ~isempty(comp_idx4)
%         plot(comp_freq.freq,squeeze(mean(20.*log10(comp_freq.powspctrm(:,comp_idx4,:)),1)),'k');
%     end;
%     axis tight;
%     xlabel('Frequency [Hz]');
%     ylabel('Power [log10]');
%     title('Power Spectrum of ICs');
%     saveas(gcf,['/data/home1/froux/DATASETS/DevelopmentData/figures/ICA_cleaning/',files1(ot).name(1:5),'_comp_idx4_2.fig'],'fig');
%     close;
%     %%
%     if ~isempty(comp_idx4)
%         figure;
%         for jt = 1:length(comp_idx4)
%             subplot(ceil(length(comp_idx4)/4),4,jt);
%             hold on;
%             imagesc(SD_EP.time,1:size(SD_EP.trial,1),squeeze(SD_EP.trial(:,jt,:)));
%             plot([0 0],[1 size(SD_EP.trial,1)],'w--');
%             plot([.4 .4],[1 size(SD_EP.trial,1)],'w--');
%             plot([1.6 1.6],[1 size(SD_EP.trial,1)],'w--');
%             axis tight;
%             xlabel('Time [s]');
%             ylabel('Trial #');
%             title(['IC # :',num2str(comp_idx4(jt))]);
%             cb = colorbar;
%             zlab = get(cb,'YLabel');
%             set(zlab,'String','IC-Amplitude [a.u]');
%         end;
%         saveas(gcf,['/data/home1/froux/DATASETS/DevelopmentData/figures/ICA_cleaning/',files1(ot).name(1:5),'_comp_idx4_3.fig'],'fig');
%         close;
%         cfg = [];
%         cfg.parameter = 'avg';
%         cfg.layout = 'CTF275.lay';
%         cfg.comment = 'no';
%         cfg.electrodes = 'off';
%         dummy = struct;
%         dummy.dimord = 'chan_time';
%         dummy.label = meg_data.label;
%         dummy.time = 0;
%
%         figure;
%         set(gcf,'Position',[50 130 1730 820]);
%         k = 0;
%         for jt = 1:length(comp_idx4)
%             k = k+1;
%             subplot(ceil(length(comp_idx4)/5),5,k);
%             dummy.avg = comp.topo(:,comp_idx4(jt));
%             ft_topoplotER(cfg,dummy);
%             title(['component #',num2str(comp_idx4(jt))]);
%         end;
%         saveas(gcf,['/data/home1/froux/DATASETS/DevelopmentData/figures/ICA_cleaning/',files1(ot).name(1:5),'_comp_idx4_4.fig'],'fig');
%         close;
%     end;