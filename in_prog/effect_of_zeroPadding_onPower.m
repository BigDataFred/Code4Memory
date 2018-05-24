%%
Fs = 1000;

t = 0:1/Fs:1;
sig1 = sin(2*pi*5.*t);

t = 0:1/Fs:10;
sig2 = sin(2*pi*5.*t);

dt = length(sig2)-length(sig1);
zp = zeros(1,dt/2);

sig1 = [zp sig1 zp];

nfft = 2^nextpow2(length(sig1));
y1 = fft(sig1,nfft)./Fs;

nfft = 2^nextpow2(length(sig2));
y2 = fft(sig2,nfft)./Fs;

y1 = y1.*conj(y1);
y2 = y2.*conj(y2);
y1 = fftshift(y1);
y2 = fftshift(y2);
y1 = y1(nfft/2:end);
y2 = y2(nfft/2:end);

f = Fs/2*linspace(0,1,nfft/2+1);

figure;
subplot(121);
plot(f,y1);
axis tight;
xlim([0 20]);
title(length(y1));
legend('data:1s');

subplot(122);
plot(f,y2,'r');
axis tight;
xlim([0 20]);
title(length(y2));
legend('data:10s');