%%
figure;
subplot(421);
hold on;
plot(spectrumLFP1{1}.fx1,10*log10(mean(spectrumLFP1{1}.S1,2)),'b','LineWidth',3);
plot(spectrumLFP2{1}.fx1,10*log10(mean(spectrumLFP2{1}.S1,2)),'r','LineWidth',3);
axis tight;

subplot(422);
hold on;
plot(freqBinz1,cleanSpectrum1,'LineWidth',3);
plot(freqBinz2,cleanSpectrum2,'r','LineWidth',3);
axis tight;

subplot(423);
hold on;
plot(spectrumSTA1{1}.fx1,10*log10(mean(spectrumSTA1{1}.S1,2)),'LineWidth',3);
plot(spectrumSTA2{1}.fx1,10*log10(mean(spectrumSTA2{1}.S1,2)),'r','LineWidth',3);

subplot(424);
hold on;
plot(freqBinz3,cleanSpectrum3,'LineWidth',3);
plot(freqBinz4,cleanSpectrum4,'r','LineWidth',3);

subplot(425);
hold on;
plot(spectrumSFC1{1}.fx1,mean(spectrumSFC1{1}.SFC1,2),'LineWidth',3);
%plot(spectrumSFC1{1}.fx1,spectrumSFC1{1}.CErr1,'r','LineWidth',3);

plot(spectrumSFC2{1}.fx1,mean(spectrumSFC2{1}.SFC1,2),'r','LineWidth',3);
axis tight;

subplot(427);
hold on;
plot(t,mean(s1,1),'LineWidth',3);
%plot(t,2*mean(E1)*ones(1,length(t)),'k');
%plot(t,-2*mean(E1)*ones(1,length(t)),'k');

plot(t,mean(s2,1),'r','LineWidth',3);

subplot(426);
hold on;
plot(spectrumSTA1{1}.fx1,((mean(spectrumSTA1{1}.S1,2))./mean(spectrumLFP1{1}.S1,2)).*100,'LineWidth',3);
plot(spectrumSTA2{1}.fx1,((mean(spectrumSTA2{1}.S1,2))./mean(spectrumLFP2{1}.S1,2)).*100,'r','LineWidth',3);

%% 
M1 = mean(SFC1,1);
M2 = mean(SFC2,1);
SE1 = std(SFC1,0,1)./sqrt(size(SFC1,1)-1);
SE2 = std(SFC2,0,1)./sqrt(size(SFC2,1)-1);
figure;
hold on;
plot(spectrumSTA1{1}{1}.fx1,M1+SE1,'b');
plot(spectrumSTA1{1}{1}.fx1,M1-SE1,'b');
plot(spectrumSTA1{1}{1}.fx1,M2+SE2,'r');
plot(spectrumSTA1{1}{1}.fx1,M2-SE2,'r');

%%
foi = 2.^([6:2:64]./6);

figure;
for it = 1:length(pvalPL)    
    
    subplot(5,5,it);
    hold on;
    plot(log10(foi),log(pvalPL{it}),'kx-');
    plot(log10([foi(1) foi(end)]),log([.05 .05]./6),'r-');
    axis tight;
    set(gca,'XTick',log10([2 3 7 13 27 53 107]));
    set(gca,'XTickLabel',[2 3 7 13 27 53 107]);
end;

%%
c = [];
for it = 1:length(spectralData.Cgrm)
    %subplot(length(spectralData.Cgrm),1,it);
    figure;
    a(it) = gca;
    hold on;
    imagesc(spectralData.tx-1,spectralData.fx,spectralData.Cgrm{it}');axis xy;
    c(it,:) = caxis;
    plot([0 0],[min(spectralData.fx) max(spectralData.fx)],'Color',[.75 .75 .75]);
    plot([2 2],[min(spectralData.fx) max(spectralData.fx)],'Color',[.75 .75 .75]);
    axis tight;
end;
caxis([0 max(max(c))]);

%%
c = [];
for it = 1:length(spectralData.SgmLFP)
    
    figure;
    %subplot(length(spectralData.SgmLFP),1,it);
    a(it) = gca;
    hold on;
    imagesc(spectralData.tx-1,spectralData.fx,spectralData.SgmLFP{it}');axis xy;
    c(it,:) = caxis;
    plot([0 0],[min(spectralData.fx) max(spectralData.fx)],'Color',[.75 .75 .75]);
    plot([2 2],[min(spectralData.fx) max(spectralData.fx)],'Color',[.75 .75 .75]);
    axis tight;
end;
%caxis([0 max(max(c))]);
%%

for it = 1:length(spectralData.SgmSPK)
    %subplot(length(spectralData.SgmSPK),1,it);
    figure;
    a(it) = gca;
    hold on;
    imagesc(spectralData.tx-1,spectralData.fx,spectralData.SgmSPK{it}'./spectralData.R{it}(:,1:31)');axis xy;
    c(it,:) = caxis;
    plot([0 0],[min(spectralData.fx) max(spectralData.fx)],'Color',[.75 .75 .75]);
    plot([2 2],[min(spectralData.fx) max(spectralData.fx)],'Color',[.75 .75 .75]);
    axis tight;
end;
%caxis([0 max(max(c))]);

%%
figure;
for it = 1:length(hitsSTA)
   subplot(2,6,it);
   hold on;
   plot(-.499:1e-3:.499,mean(hitsSTA{it},1),'k');
   plot(-.499:1e-3:.499,mean(missSTA{it},1),'g');
end;

%%
figure;
k = 0;
for it = 1:size(LFP2FLPcoh.C,1)    
    for jt = 1:size(LFP2FLPcoh.C,2)
        k = k+1;
        subplot(4,3,k);
        plot(LFP2FLPcoh.fx,squeeze(LFP2FLPcoh.C{it,jt}));
        axis xy;
    end;
end;