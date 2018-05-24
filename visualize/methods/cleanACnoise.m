
%%
Fs = 1000;
t = 0:1/Fs:7;

sig = sin(2*pi*10.*t)+sin(2*pi*50.*t)+0+0.1*randn(1,length(t));

[cleanSignal,~] = CleanLineNoise( sig ,'Fs', Fs , 'noiseFreq', 50,'windowSize',1);

figure;
subplot(221);
hold on;
plot(t,sig);
plot(t,cleanSignal,'r');
xlim([0 .1]);

subplot(222);
nfft = 2^nextpow2(length(sig));
y = fft(sig,nfft)./Fs;
y = fftshift(y);
y = y.*conj(y);
y = y(nfft/2:end);

y2 = fft(cleanSignal,nfft)./Fs;
y2 = fftshift(y2);
y2 = y2.*conj(y2);
y2 = y2(nfft/2:end);

f = Fs/2*linspace(0,1,nfft/2+1);

hold on;
plot(f,y);
plot(f,y2,'r');
xlim([0 100]);

subplot(223);
plot(t,sin(2*pi*10.*t)+cleanSignal+sin(2*pi*50.*t),'r');
xlim([0 .1]);