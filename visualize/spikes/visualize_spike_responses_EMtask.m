%%
figure;
subplot(321);
hold on;
plot([0 0],[1 1],'r^','MarkerFaceColor','r');
h = area([0 2e3],[length(trl_all)+1 length(trl_all)+1],0);
set(h,'FaceColor','r','FaceAlpha',.2,'EdgeColor','r');

plot([2e3 2e3],[1 1],'bv','MarkerFaceColor','b');
h = area([2e3 4e3],[length(trl_all)+1 length(trl_all)+1],0);
set(h,'FaceColor','b','FaceAlpha',.2,'EdgeColor','b');

for kt = 1:length(ts_trl)
    x = ts_trl{kt};
    x = x(x>=dt(1) & x <=dt(length(dt)));
    y = kt*ones(1,length(x));
    x = [x;x];
    y = [y-.5;y+.5];
    line(x,y,'Color','k');
end;
axis tight;xlim([dt(1) dt(end)]);
title(length([ts_trl{:}]));

subplot(323);
hold on;
%Y = sum(n2,2)/size(n2,2)/((dt3(2)-dt3(1))/Fs);
Y = mean(n2,2)/((dt3(2)-dt3(1))/Fs);
SE = (std(n2,0,2)./sqrt(size(n2,2)-1))./((dt3(2)-dt3(1))/Fs);
plot([0 0],[min(Y-SE) max(Y+SE)],'r');
plot([2e3 2e3],[min(Y-SE) max(Y+SE)],'b');
plot(dt3,Y,'-','Color',[.5 .5 .5]);
errorbar(dt3,Y,SE,'s','Color',[.5 .5 .5]);
plot([0 0],[min(Y-SE) min(Y-SE)],'r^','MarkerFaceColor','r');
plot([2e3 2e3],[min(Y-SE) min(Y-SE)],'bv','MarkerFaceColor','b');
axis tight;
xlim([dt(1) dt(end)]);
title(round(mean(mean(n))/1e-3*1e2)/1e2);

subplot(322);
bar(dt2,n_isi);
axis tight;
title(pct);

%                     subplot(524);
%                     hold on;
%                     plot(dt2,cumsum(n_isi_b./sum(n_isi_b)),'k');
%                     plot(dt2,cumsum(n_isi_c./sum(n_isi_c)),'r');
%                     plot(dt2,cumsum(n_isi_e./sum(n_isi_e)),'b');
%                     axis tight;ylim([0 1]);

subplot(324);
hold on;
plot(spike_dat.waveformtime,wvf,'r');
plot(spike_dat.waveformtime,mean(wvf,1),'k','LineWidth',3);
axis tight;


figure;
subplot(321);
hold on;
plot(-win:win,mean(STA1,2),'Color',[0 0 0]);
plot(-win:win,mean(STA2,2),'Color',[.9 0 0]);
plot(-win:win,mean(STA3,2),'Color',[0 0 .9]);
axis tight;

subplot(323);
hold on;
plot(f1,Snn1./R1,'k');
plot(f1,Snn2./R2,'r');
plot(f1,Snn3./R3,'b');
axis tight;
title(num2str((R1+R2+R3)./3));

subplot(324);
hold on;
plot(f1,(Syy1),'k');
plot(f1,(Syy2),'r');
plot(f1,(Syy3),'b');
axis tight;

subplot(325);
hold on;
plot(f2,C1,'Color','k');
plot(f2,C2,'Color','r');
plot(f2,C3,'Color','b');
axis tight;

subplot(326);
hold on;
imagesc((tAx-2).*1e3,fAx,Cgrm');
plot([0 0],[fAx(1) fAx(end)],'Color',[.85 .85 .85]);
plot([0 0],[fAx(1) fAx(1)],'r^','MarkerFaceColor','r');
plot([2e3 2e3],[fAx(1) fAx(end)],'Color',[.85 .85 .85]);
plot([2e3 2e3],[fAx(1) fAx(1)],'bv','MarkerFaceColor','b');
axis xy;axis tight;