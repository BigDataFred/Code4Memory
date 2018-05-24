%%
Fs = 1e3;
t = 0:1/Fs:1;

theta1 = -.85+sin(2*pi*5.*t)+.05*randn(1,length(t));

Fs = 1/1e-3;
TW = 1;%(size(data1,1)/Fs)*1/(size(data1,1)/Fs);
params1                  = [];
params1.tapers           = [TW , 2*TW-1];
params1.pad              = 1;
params1.Fs               = Fs;
params1.err              = 0;
params1.trialave         = 0;
params1.fpass            = [0 300];

[S1,f] = mtspectrumc(theta1',params1);

Fs = 1e3;
t = 0:1/Fs:1;

theta2 = 0+sin(2*pi*5.*t)+.05*randn(1,length(t));

Fs = 1/1e-3;
TW = 1;%(size(data1,1)/Fs)*1/(size(data1,1)/Fs);
params1                  = [];
params1.tapers           = [TW , 2*TW-1];
params1.pad              = 1;
params1.Fs               = Fs;
params1.err              = 0;
params1.trialave         = 0;
params1.fpass            = [0 300];

[S2,f] = mtspectrumc(theta2',params1);

Fs = 1e3;
t = 0:1/Fs:1;

theta3 = .85+sin(2*pi*5.*t)+.05*randn(1,length(t));

Fs = 1/1e-3;
TW = 1;%(size(data1,1)/Fs)*1/(size(data1,1)/Fs);
params1                  = [];
params1.tapers           = [TW , 2*TW-1];
params1.pad              = 1;
params1.Fs               = Fs;
params1.err              = 0;
params1.trialave         = 0;
params1.fpass            = [0 300];

[S3,f] = mtspectrumc(theta3',params1);

figure;
subplot(321);
hold on;
plot(t,theta1);
plot([t(1) t(end)],[0 0],'r');
axis tight;
subplot(322);
plot(f,S1);
xlim([0 30]);
axis tight;
xlim([0 30]);
subplot(323);
hold on;
plot(t,theta2);
plot([t(1) t(end)],[0 0],'r');
axis tight;
subplot(324);
plot(f,S2);
xlim([0 30]);
axis tight;
xlim([0 30]);
subplot(325);
hold on;
plot(t,theta3);
plot([t(1) t(end)],[0 0],'r');
axis tight;
subplot(326);
plot(f,S3);
xlim([0 30]);
axis tight;
xlim([0 30]);