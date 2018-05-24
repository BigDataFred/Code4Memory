%% simulate 1/f 
%zeroes and poles from www.firstpr.com.au/dsp/pink-noise
z = [0.98443604;...
     0.83392334;...
     0.07568359];

p = [0.99572754;...
     0.94790649;...
     0.53567505];
 
k = 0.3;
     
[b,a] = zp2tf(z,p,k);

x = randn(1,1e5);

y = filtfilt(b,a,[x x x]);
y = y(length(x)+1:length(x)*2);

Fs = 1000;
t = 0:1/Fs:(length(x)-1)/Fs;

figure;
subplot(2,2,1:2);
hold on;
plot(t,x);
plot(t,y,'r');
ylabel('a.u.');
xlabel('Time (s)');
title('Time domain');

subplot(223);
nfft = 2^nextpow2(length(x));
p1 = fft(x,nfft)./nfft;
p2 = fft(y,nfft)./nfft;
p1 = fftshift(p1);
p2 = fftshift(p2);
p1 =  p1.*conj(p1);
p2 =  p2.*conj(p2);
p1 = p1(nfft/2:end);
p2 = p2(nfft/2:end);

n = length(p1);
b = regress(log10(p2)',[ones(length(p2),1) log10(1:n)']);

f = Fs/2*linspace(0,1,nfft/2+1);


hold on;
[ax,h1,h2] = plotyy(f,p1,f,p2);
set(ax,'XLim',[0 10]);
ylabel('a.u.');
xlabel('Freq (Hz)');
title('Frequency domain');
set(ax(2),'YTick',[]);

subplot(224);
[ax,h1,h2] = plotyy(log10(f),log10(p1),log10(f),log10(p2));
set(ax,'XLim',log10([0 100]));
ylabel('log');
xlabel('log');
title('Frequency domain');
set(ax(2),'YTick',[]);

%% relationship between derivative in time domain and correcting for 1/f in frequency domain

tm = 1;% duration of data in s
Fs = 1000; % set sampling rate.
t = 0:1/Fs:tm;% time vector
w = 2*pi;% freq in radians

fx  = sin(w.*t) + sin(2*w.*t) + sin(4*w.*t) + sin(6*w.*t);
fx2 = cos(w.*t) + cos(2*w.*t) + cos(4*w.*t) + cos(6*w.*t);


nfft= 2^nextpow2(length(fx));

f = Fs/2*linspace(0,1,nfft/2+1).*w;

y = fft(fx,nfft)./Fs;
y = fftshift(y);
y = y(nfft/2:end);
y = y.*f;
y = y.*conj(y);
[y1] = y;

y = fft(fx2,nfft)./Fs;
y = fftshift(y);
y = y(nfft/2:end);
y = y.*f;
y = y.*conj(y);
[y2] = y; 

y = fft(diff(fx)./(1/Fs),nfft)./Fs;
y = y.*conj(y);
y = fftshift(y);
[y3] = y(nfft/2:end);

figure;

subplot(221);
plot(t,fx,'b');
axis tight;
legend('sin(\omegat)+sin(2*\omegat)+sin(4*\omegat)+sin(6*2\omegat)');
xlabel('t');
title('Time domain');

subplot(223);
%[ax,h1,h2] = plotyy(t(2:end),diff(fx)./(t(2)-t(1)),t,fx2);
[ax,h1,h2] = plotyy(t,gradient(fx),t,fx2);
set(h1,'Color','k','LineWidth',3);
set(h2,'Color','r');
axis tight;
legend([h1 h2],'dy/dx','cos(\omega*t)');
xlabel('t');
title('Time domain');

% subplot(222);
% hold on;
% plot(f,y1,'b');
% plot(f,y2,'r');
% xlim([0 60]);
% xlabel('w (rad)');
% legend('sin(w*t)','cos(w*t)');
% %set(gca,'XTick',[1 2 4 6]);
% title('f(w)');

subplot(224);
[ax,h1,h2] = plotyy(f,y3,f,y2);
set(h1,'Color','k','LineWidth',3);
set(h2,'Color','r');
hold(ax(2),'on');
%h3 = plot(ax(2),f,y1,'b');
set(ax,'Xlim',[0 60]);
set(ax(1),'YColor','k');
set(ax(2),'YColor','b');
xlabel('\omega');
title('Frequency domain');

legend([h1,h2],'dy/dx','cos(w*t)');
%% masking effect of 1/f
x = 1:.1:30;
x = x';
a1 = 1;
a2 = 1.1;
fx = 1./x.^a1;
fx2 = 1./x.^a2;

fx2(find(x==4.8)) = fx2(find(x==4.8))+ 0.75e-5;
fx2(find(x==5)) = fx2(find(x==5))+ 1e-3;
fx2(find(x==5.1)) = fx2(find(x==5.1))+ 0.75e-5;

figure;
subplot(221);
hold on;
plot(x,fx);
plot(x,fx2,'r');
legend('baseline:1/f^{a1}','poststim:1/f^{a2}');
xlabel('Freq. (Hz)');
ylabel('Power (a.u.)');
axis tight;
title('See an oscillation?');

subplot(222);
hold on;
plot(log10(x),log10(fx));
plot(log10(x),log10(fx2),'r');

b1 = regress(log10(fx),[ones(length(fx),1) log10(x)]);
yp1 = b1(1)+b1(2)*log10(x);
c1 = log10(fx)-yp1;

b2 = regress(log10(fx2),[ones(length(fx),1) log10(x)]);
yp2 = b2(1)+b2(2)*log10(x);
c2 = log10(fx2)-yp2;
legend(['a1:',num2str(round(b1(2)*10)/10)],['a2:',num2str(round(b2(2)*10)/10)]);
xlabel('Freq. (log)');
ylabel('Power (log)');
axis tight;

subplot(223);
hold on;
plot(x,10.^c1);
plot(x,10.^c2,'r');
legend('baseline','poststim');
xlabel('Freq. (Hz)');
ylabel('Power (a.u.)');
axis tight;
title('How about now?');

