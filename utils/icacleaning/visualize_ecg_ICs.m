function visualize_ecg_ICs(meg_data,comp,comp_idx1,qrsA,comp_freq)

if nargin == 0;
    comp_idx1 = [];
    qrsA = {};
end;       
%%
cfg = [];
cfg.layout = 'neuromag306planar.lay';
cfg.comment = 'no';
cfg.electrodes = 'off';

dum = struct;
dum.time = 0;
dum.dimord = 'chan_time';
dum.label = meg_data.label;
%% visualize results of automatic EKG detection
if ~isempty(comp_idx1)               
    k = 0;
    figure;
    set(gcf,'Position',[50 250 1730 620]);
    for jt = 1:length(comp_idx1)
        
        k = k+1;
        subplot(length(comp_idx1),5,k);
        dum.avg = comp.topo(:,comp_idx1(jt));
        ft_topoplotER(cfg,dum);
        title(['component #',num2str(comp_idx1(jt))]);
        
        k = k+1;
        subplot(length(comp_idx1),5,k);
        hold on;
        plot(comp_freq.freq,squeeze(log10(comp_freq.powspctrm(:,comp_idx1(jt),:))),'Color',[.75 .75 .75]);
        plot(comp_freq.freq,squeeze(mean(log10(comp_freq.powspctrm(:,comp_idx1(jt),:)),1)),'k');
        axis tight;
        xlabel('Frequency [Hz]');
        ylabel('Power [log10]');
        title('Power Spectrum');
        
        k = k+1;
        subplot(length(comp_idx1),5,k);
        hold on;
        dums = comp.trial{1}(comp_idx1(jt),:);
        dums = dums-mean(dums);
        dums = dums./std(dums);
        plot(comp.time{1},dums);
        plot([min(comp.time{1}) max(comp.time{1})],[2 2],'r--');
        plot([min(comp.time{1}) max(comp.time{1})],[-2 -2],'r--');
        axis tight;
        xlabel('Time [s]');
        ylabel('Amplitude [\sigma]');
        title('Time course');        
        
        k = k+1;
        subplot(length(comp_idx1),5,k);
        hold on;
        [ax,h1,h2] = plotyy((0:size(qrsA{jt}.avg,2)-1)./comp.fsample,qrsA{jt}.avg,(0:size(qrsA{jt}.avg,2)-1)./comp.fsample,qrsA{jt}.itc);
        axis(ax,'tight');
        xlabel(ax(1),'Time [s]');
        ylabel(ax(1),'Power [\sigma]');
        ylabel(ax(2),'ITC [\sigma]');
        title('QRS-locked average');

        k = k+1;
        subplot(length(comp_idx1),5,k);
        imagesc((0:size(qrsA{jt}.avg,2)-1)./comp.fsample,1:size(qrsA{jt}.trl,1),(abs(qrsA{jt}.trl)));axis xy;
        cb = colorbar;
        zlab = get(cb,'YLabel');
        set(zlab,'String','[a.u.]');
        axis tight;
        xlabel('Time [s]');
        ylabel('Trial  #');
        title('Temporal distribution across trials');
    end;    
end;