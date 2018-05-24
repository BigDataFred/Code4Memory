%%
Fs = 1000;
t = 0:1/Fs:1;
sig = sin(2*pi*1.*t)+sin(2*pi*5.*t)+sin(2*pi*10.*t);%+(0+1*randn(1,length(t)));

sig2 = [zeros(1,length(t)*5) sig zeros(1,length(t)*5)];

nfft = 2^nextpow2(length(sig));
y = fft(sig,nfft)./Fs;
y = fftshift(y);
y = y(nfft/2:end);
y = y.*conj(y);
f = Fs/2*linspace(0,1,nfft/2+1);

nfft2 = 2^nextpow2(length(sig2));
y2 = fft(sig,nfft2)./Fs;
y2 = fftshift(y2);
y2 = y2(nfft2/2:end);
y2 = y2.*conj(y2);
f2 = Fs/2*linspace(0,1,nfft2/2+1);

figure;
hold on;
plot(f,y,'ro-','LineWidth',2);
plot(f2,y2,'b.','LineWidth',2);

xlim([0 50]);