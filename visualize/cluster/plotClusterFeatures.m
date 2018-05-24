%%
figure;
plot(spike_dat.waveformtime,wvf);
title(['SNR : ',num2str(snr)]);

% figure;
% plot(spike_dat.waveformtime,itc);
% 
% figure;
% plot(lag,Cxx);
% title(['Percent :',num2str(pctZero)]);
% 
% figure;
% hold on;
% plot(spike_stats.f,spike_stats.Pxxn);
% %plot([spike_stats.f(1) spike_stats.f(end)],[m+10*s m+10*s],'r--');
% 
% dt = 0:1:1000;
% [n,~] = hist(isi,dt);
% 
% parm1 = gamfit(isi);
% yGam = gampdf(dt, parm1(1),parm1(2));
% yGam = yGam*length(isi);
% 
% parm2 = poissfit(isi);
% yPois = poisspdf(dt, parm2);
% yPois = yPois*length(isi);
% 
% figure;
% hold on;
% bar(dt,n);
% plot(dt,yGam,'r');
% plot(dt,yPois,'c');
% title(['Percent < 3ms:',num2str(pctBelow),'Percent > 1000ms:',num2str(pctAbove)]);

