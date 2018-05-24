
%%
Fs1 = 32e3;
t = 0:1/Fs1:100;
[sig1] = sin(2*pi*50.*t) + .15*randn(1,length(t));

[cleanSignal1,~] = CleanLineNoise(sig1 ,'Fs', Fs1 , 'noiseFreq', 50,'windowSize',1);
[cleanSignal2,~] = CleanLineNoise(sig1 ,'Fs', Fs1 , 'noiseFreq', 50,'windowSize',90);

Fs2 = 1e3;
t = 0:1/Fs2:100;
[sig2] = sin(2*pi*50.*t) + .15*randn(1,length(t));

[cleanSignal3,~] = CleanLineNoise(sig2 ,'Fs', Fs2 , 'noiseFreq', 50,'windowSize',1);
[cleanSigna4,~] = CleanLineNoise(sig2 ,'Fs', Fs2 , 'noiseFreq', 50,'windowSize',100);


params                  = [];
params.pad              = 2;
params.fpass            = [0 100];
params.tapers           = [3 5];
params.Fs               = Fs1;

[S1,f1] = mtspectrumc(sig1',params);
[S2,f1] = mtspectrumc(cleanSignal1',params);
[S3,f1] = mtspectrumc(cleanSignal2',params);

params.Fs               = Fs2;

[S4,f2] = mtspectrumc(sig2',params);
[S5,f2] = mtspectrumc(cleanSignal3',params);
[S6,f2] = mtspectrumc(cleanSigna4',params);

%%
figure;
hold on;
plot(f1,S1,'b');
plot(f1,S2,'r');
plot(f1,S3,'k');

figure;
hold on;
plot(f2,S4,'b');
plot(f2,S5,'r');
plot(f2,S6,'k');

%%
ix = find(hPL==1);

figure;
for it = 1:length(spectralData.SgmLFP)
    subplot(2,4,it);
    imagesc(spectralData.tx,spectralData.fx,squeeze(spectralData.Cgrm{it})');
    axis xy;
    title(chID(ix(it)));
end;
