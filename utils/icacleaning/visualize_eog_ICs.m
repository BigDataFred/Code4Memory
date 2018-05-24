function visualize_eog_ICs(eog_idx,comp,comp_freq,bData)
%%
cfg = [];
cfg.layout = 'neuromag306planar.lay';
cfg.comment = 'no';
cfg.electrodes = 'off';

k = 0;
figure;
set(gcf,'Position',[50 250 1730 620]);
X = zeros(length(eog_idx)*length(comp.trial),length(comp.time{1}));
for jt = 1:length(eog_idx)
    k = k+1;
    subplot(length(eog_idx)+1,4,k);
    cfg.component = eog_idx(jt);
    ft_topoplotIC(cfg,comp);
    title(['component #',num2str(eog_idx(jt))]);
    
    k = k+1;
    subplot(length(eog_idx)+1,4,k);
    hold on;
    plot(comp_freq.freq,squeeze(log10(comp_freq.powspctrm(:,eog_idx(jt),:))),'Color',[.75 .75 .75]);
    h = [];
    h = plot(comp_freq.freq,squeeze(mean(log10(comp_freq.powspctrm(:,eog_idx(jt),:)),1)),'k');
    axis tight;
    xlabel('Frequency [Hz]');
    ylabel('Power [log10]');
    
    k = k+1;
    subplot(length(eog_idx)+1,4,k);
    hold on;
    E = zeros(1,length(comp.trial));
    for kt = 1:length(comp.trial)
        E(kt) = sqrt(mean(abs(hilbert(comp.trial{kt}(eog_idx(jt),:))).^2));
    end;
    E = abs((E-mean(E))./std(E));
    
    dums = comp.trial{find(E==max(E))}(eog_idx(jt),:);
    dums = dums-mean(dums);
    dums = dums./std(dums);
    plot(comp.time{1},dums);
    plot([min(comp.time{1}) max(comp.time{1})],[2 2],'r--');
    plot([min(comp.time{1}) max(comp.time{1})],[-2 -2],'r--');
    axis tight;
    xlabel('Time [s]');
    ylabel('IC-Amplitude [\sigma]');
    
    k = k+1;
    subplot(length(eog_idx)+1,4,k);
    trlmat = zeros(length(comp.trial),length(comp.time{1}));
    for kt = 1:length(comp.trial)
        trlmat(kt,:) = comp.trial{kt}(eog_idx(jt),:);
        trlmat(kt,:) = (trlmat(kt,:)-mean(trlmat(kt,:)))./std(trlmat(kt,:));
    end;
    
    imagesc(comp.time{1},1:length(comp.trial),abs(trlmat)>4);
    caxis([-1 1]);
    axis tight;
    xlabel('Time [s]');
    ylabel('Trial  #');
    
    X = [X;trlmat];
    
end;

k = k+1;
% subplot(length(eog_idx)+1,4,k);
% hold on;
% plot(eogcomp_coh.freq,eogcomp_coh.cohspctrm);
% plot([min(eogcomp_coh.freq) max(eogcomp_coh.freq)],[4 4],'r--');
% axis tight;
% xlabel('Frequency [Hz]');
% ylabel('Coherence [\sigma]');
% title('EOG-ICs coherence spectrum');

k = k+1;
subplot(length(eog_idx)+1,4,k);
hold on;

for jt = 1:length(bData.bcMEG)
    x2 = squeeze(mean(bData.bcMEG{jt}.trl(eog_idx,:,:),2));
    for kt = 1:size(x2,1)
        x2(kt,:) = x2(kt,:)-min(x2(kt,:));
        x2(kt,:) = x2(kt,:)./(max(x2(kt,:))-min(x2(kt,:)));
    end;
end;
h = [];
h = plot(linspace(-.5,.5,length(x2)),x2);
xlabel('Time [s]');
ylabel('[a.u.]');
legend(h,num2str(eog_idx));


% subplot(length(eog_idx)+1,4,k);
% 
% plot(filt_eeg_data1.time{1},filt_eeg_data1.trial{find(E==max(E))});
% axis tight;
% xlabel('Time [s]');
% ylabel('EOG-Amplitude [\muV]');

k = k+1;
subplot(length(eog_idx)+1,4,k);
hold on;
br = conv(sum(abs(X)>4,1),hanning(40));
br = br./sum(hanning(40));
br = br(20:end-20);
br = (br-mean(br))./std(br);
plot(comp.time{1},br,'k.');
h = [];
h = plot([min(comp.time{1}) max(comp.time{1})],[2 2],'r--');
axis tight;
xlabel('Time [s]');
ylabel('Count [\sigma]');
legend(h,'\sigma = 2');

