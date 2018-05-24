%%
load P02_fVSpEM_2016-07-15_11-16-26_lfpDataStimLockedSegmentedAVGdownsampled.mat;

it = 1;

[nsamp,ntrl] = size(LFPavg{1});

T = nsamp/Fs;
W = 3/T;
TW = T*W;
k = 2*TW-1;
params                  = [];
params.Fs               = Fs;
params.pad              = 8;
params.fpass            = [0 200];
params.trialave         = 0;
params.tapers           = [TW k];

[S,f] = mtspectrumc( (LFPavg{it}), params );

ix = find(f>=45 & f <=55);
for jt = 1:ntrl   
    fprintf([num2str(jt),'/',num2str(ntrl)]);
    [~,m] = max(S(ix,jt));
    ix = ix(m);
    
    [LFPavg{it}(:,jt),noise] = CleanLineNoise(LFPavg{it}(:,jt)','Fs',Fs,'noiseFreq',f(ix),'windowSize',nsamp/Fs);
    fprintf('\n');
end;

[S2,f] = mtspectrumc( (LFPavg{it}), params );

figure;
hold on;
plot(f,mean(log10(S),2));
plot(f,mean(log10(S2),2),'r');

%%
[dum,noise] = CleanLineNoise(x,'Fs',Fs,'noiseFreq',f(ix),'windowSize',3.5);

T = length(x)/Fs;
W = 1/T;
TW = T*W;
k = 2*TW-1;
params                  = [];
params.Fs               = Fs;
params.pad              = 0;
params.fpass            = [0 100];
params.trialave         = 0;
params.tapers           = [TW k];

[S2,f] = mtspectrumc( dum', params );
[S3,f] = mtspectrumc( noise', params );

figure;
subplot(221);
hold on;
plot(x(1:3/50*Fs));
plot(dum(1:3/50*Fs),'r');
axis tight;
title('LFP signal');
legend('raw LFP','cleaned LFP');
xlabel('Sample #');
ylabel('Amplitude [\muV]');
subplot(222);
plot(noise(1:3/50*Fs));
axis tight;
title('fitted noise signal');
xlabel('Sample #');
ylabel('Amplitude [\muV]');
subplot(223);
hold on;
plot(f,20*log10(S));
plot(f,20*log10(S2),'r');
legend('raw LFP','cleaned LFP');
axis tight;
xlabel('Frequency [Hz]');
ylabel('Power [dB]');
title('LFP signal');
subplot(224);
hold on;
plot(f,20*log10(S3));
axis tight;
xlabel('Frequency [Hz]');
ylabel('Power [dB]');
title('fitted noise signal');
