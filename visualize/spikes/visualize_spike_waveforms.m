function visualize_spike_waveforms(spike)
tx = -5e-4:(1/spike.hdr.SamplingFrequency):4.99e-4;
tx =tx.*1e2;

scrsz = get(groot,'ScreenSize');
figure('Position',[1 0 scrsz(3)/3 scrsz(4)]);
subplot(4,1,1:3);
a = gca;
imagesc(tx,1:length(spike.unit{1}),squeeze(spike.waveform{1}(1,:,:))');
axis xy;
subplot(4,1,4);
b = gca;
hold on;
plot(tx,squeeze(spike.waveform{1}(1,:,:)),'Color',[.75 .75 .75]);
plot(tx,mean(squeeze(spike.waveform{1}(1,:,:)),2),'r','LineWidth',3);
axis tight;

xlabel(a,'Time [ms]');
xlabel(b,'Time [ms]');

ylabel(a,'Number of detected spikes');
ylabel(b,'Amplitude [\muV]');