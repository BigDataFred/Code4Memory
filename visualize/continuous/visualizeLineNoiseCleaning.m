%%
s = 5;
[step] = round(length(dataSamples)/(s*Fs));

ix = 1:Fs*s;
dum = zeros(1,length(dataSamples));

for it = 1:step
    
    fprintf([num2str(it),'/',num2str(step)]);
    
    T = length(dataSamples(ix))/Fs;
    W = 1/T;
    TW = round(T*W);
    k = round(2*TW-1);
    params                  = [];
    params.Fs               = Fs;
    params.pad              = 6;
    params.fpass            = [40 60];
    params.trialave         = 0;
    params.tapers           = [TW k];
    
    [S,f] = mtspectrumc( dataSamples(ix)', params );
    
    ix2 = find(f>=45 & f <=55);
    [~,m] = max(S(ix2));
    ix2 = ix2(m);
    
    [dum(ix),noise] = CleanLineNoise(dataSamples(ix),'Fs',Fs,'noiseFreq',f(ix2),'windowSize',s);
    
%     [Sn,fn] = mtspectrumc( noise', params );
%     
%     [b,a] = butter(2,[49.5 50.5]./(Fs/2),'stop');% apply band-stop for LFP
%     [dum2] = filtfilt(b,a,dataSamples(ix));
%     
%     params                  = [];
%     params.Fs               = Fs;
%     params.pad              = 0;
%     params.fpass            = [40 60];
%     params.trialave         = 0;
%     params.tapers           = [TW k];
%     
%     [S1,f1] = mtspectrumc( dataSamples(ix)', params );
%     [S2,f2] = mtspectrumc( dum(ix)', params );
%     [S3,f3] = mtspectrumc( dum2', params );
%     
%     params                  = [];
%     params.Fs               = Fs;
%     params.pad              = 0;
%     params.fpass            = [0 200];
%     params.trialave         = 0;
%     params.tapers           = [TW k];
%     
%     [S4,f4] = mtspectrumc( dum(ix)', params );
%     
%     figure;
%     subplot(421);
%     hold on;
%     plot(f,20*log10(S),'LineWidth',3);
%     plot(fn,20*log10(Sn),'r');
%     axis tight;
%     xlabel('Frequency [hz]');
%     ylabel('Power [dB]');
%     title('LFP and noise spectra');
%     
%     subplot(422);
%     plot([1:3/50*Fs]./Fs,noise(1:3/50*Fs));
%     axis tight;
%     title(['Line freq: ',num2str(f(ix2)),'Hz']);
%     xlabel('Time [s]');
%     ylabel('Amplitude [a.u]');
%     
%     subplot(423);
%     hold on;
%     plot(f1,20*log10(S1));
%     plot(f3,20*log10(S3),'r');
%     xlim([48 52]);
%     xlabel('Frequency [hz]');
%     ylabel('Power [dB]');
%     title('50 Hz notch filter');
%     
%     subplot(424);
%     hold on;
%     plot(f1,20*log10(S1));
%     plot(f2,20*log10(S2),'r');
%     xlim([48 52]);
%     xlabel('Frequency [hz]');
%     ylabel('Power [dB]');
%     title('Template subtraction');
%     
%     subplot(4,2,5:6);
%     hold on;
%     plot(ix./Fs- min(ix./Fs),dataSamples(ix),'LineWidth',3);
%     plot(ix./Fs- min(ix./Fs),dum(ix),'r');
%     axis tight;
%     xlim([8 8.1]);
%     xlabel('Time [s]');
%     ylabel('Amplitude [a.u]');
%     title('LFP pre/post cleaning');
%     
%     subplot(427);
%     hold on;
%     plot(f4,20*log10(S4),'r');
%     xlabel('Frequency [hz]');
%     ylabel('Power [dB]');
%     title('LFP Spectrum after cleaning');
    
    ix = ix+s*Fs;
    if max(ix) > length(dataSamples)
        ix = ix(1):length(dataSamples);
    end;
    
    fprintf('\n');
    
end;

%%
[b,a] = butter(2,[49.5 50.5]./(Fs/2),'stop');%
[dum2] = filtfilt(b,a,dataSamples);

T = length(dataSamples)/Fs;
W = 1/T;
TW = round(T*W);
k = round(2*TW-1);

params                  = [];
params.Fs               = Fs;
params.pad              = 0;
params.fpass            = [0 200];
params.trialave         = 0;
params.tapers           = [TW k];

[S,f] = mtspectrumc( dataSamples', params );
[S2,f2] = mtspectrumc( dum', params );
[S3,f3] = mtspectrumc( dum2', params );

figure;
subplot(421);
hold on;
plot(f,20*log10(S),'LineWidth',3);
plot(f3,20*log10(S3),'r');
axis tight;
xlabel('Frequency [hz]');
ylabel('Power [dB]');
title('50Hz Notch filter');
subplot(422);
hold on;
plot(f,20*log10(S),'LineWidth',3);
plot(f2,20*log10(S2),'r');
axis tight;
xlabel('Frequency [hz]');
ylabel('Power [dB]');
title('Template subtraction');
subplot(423);
hold on;
plot(f,20*log10(S),'LineWidth',3);
plot(f3,20*log10(S3),'r');
axis tight;
xlim([48 52]);
xlabel('Frequency [hz]');
ylabel('Power [dB]');
subplot(424);
hold on;
plot(f,20*log10(S),'LineWidth',3);
plot(f2,20*log10(S2),'r');
axis tight;
xlabel('Frequency [hz]');
ylabel('Power [dB]');
xlim([48 52]);
subplot(425);
hold on;
plot(f,20*log10(S),'LineWidth',3);
plot(f3,20*log10(S3),'r');
axis tight;
xlim([98 102]);
xlabel('Frequency [hz]');
ylabel('Power [dB]');
subplot(426);
hold on;
plot(f,20*log10(S),'LineWidth',3);
plot(f2,20*log10(S2),'r');
axis tight;
xlabel('Frequency [hz]');
ylabel('Power [dB]');
xlim([98 102]);
subplot(427);
hold on;
plot(f,20*log10(S),'LineWidth',3);
plot(f3,20*log10(S3),'r');
axis tight;
xlim([148 152]);
xlabel('Frequency [hz]');
ylabel('Power [dB]');
subplot(428);
hold on;
plot(f,20*log10(S),'LineWidth',3);
plot(f2,20*log10(S2),'r');
axis tight;
xlabel('Frequency [hz]');
ylabel('Power [dB]');
xlim([148 152]);

%%
selIx = selIx7;
dt = [1*Fs 5*Fs];
tx = [-dt(1):dt(2)]./Fs;

x = [];
x2 = [];
x3 = [];
for it =1:length(selIx)    
    x(it,:) = dataSamples(TTLix(selIx(it))-dt(1):TTLix(selIx(it))+dt(2));
    x2(it,:) = dum(TTLix(selIx(it))-dt(1):TTLix(selIx(it))+dt(2));
    x3(it,:) = dum2(TTLix(selIx(it))-dt(1):TTLix(selIx(it))+dt(2));
end;

movingwin               = [0.25 0.001];
T = movingwin(1);
W = 10;
TW = T*W;
k = 2*TW-1;

params                  = [];
params.pad              = 0;
params.fpass            = [30 200];
params.Fs               = Fs;
params.tapers           = [TW k];
params.trialave         = 1;

[S1,t,f] = mtspecgramc( gradient(x)', movingwin, params );
[S2,t,f] = mtspecgramc( gradient(x2)', movingwin, params );
[S3,t,f] = mtspecgramc( gradient(x3)', movingwin, params );

figure;
imagesc(t-1,f,10*log10(S1)');
axis xy;
ylabel('Frequency [Hz]');
xlim([-.4 4]);
title('LFP before line noise cleaning');

figure;
imagesc(t-1,f,10*log10(S2)');
axis xy;
ylabel('Frequency [Hz]');
xlim([-.4 4]);
title('LFP after line noise cleaning');

figure;
imagesc(t-1,f,10*log10(S3)');
axis xy;
ylabel('Frequency [Hz]');
xlim([-.4 4]);
title('LFP after line noise cleaning');
