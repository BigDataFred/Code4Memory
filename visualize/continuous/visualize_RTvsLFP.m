function visualize_RTvsLFP(tlck,LogDat)
%%
n = size(tlck.trial,1);

figure;
subplot(4,1,1:3);
hold on;
imagesc(tlck.time,1:n,squeeze(tlck.trial));
axis xy;
plot(zeros(1,n),[1:n],'k.');
plot(LogDat.RT,1:n,'m.');
axis tight;
caxis([-50 50]);

subplot(4,1,4);
plot(tlck.time,tlck.avg);
xlim([min(tlck.time) max(LogDat.RT)]);