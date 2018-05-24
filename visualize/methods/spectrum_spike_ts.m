%%
clear;
clc;
close all;

addpath(genpath('/media/rouxf/rds-share/Fred/code/mcode/toolboxes/chronux_2_11/'));

%%
Fs = 1000;
t = 1e-4:1/Fs:5;
x = zeros(1,length(t));
x2 = sin(2*pi*5.*t);

sel = find(x2 >= max(x2)/100*95);

x(sel) =1;

sel = randperm(length(t));
sel = sel(1:length(t)/2);

x(sel) =1;

%%
ts = t(find(x==1));
dt = min(t):1e-3:max(t);
n = histc(ts,dt);

[xc,lag] = xcorr(n-mean(n),'coeff');

nfft = 2^nextpow2(length(xc));
y= fft(xc,nfft)./Fs;
y = conj(y).*y;
y = fftshift(y);
y = y(nfft/2:end);

f = Fs/2*linspace(0,1,nfft/2+1);

T = t(end);
W = 1/T;

params = [];
params.tapers = [W T 1];
params.Fs = Fs;
params.pad = 0;
params.fpass = [0 100];
params.err = [2 0.05];
params.trialave = 1;

[S,f2,Serr] = mtspectrumpb(n,params,0);

figure;
subplot(221);
plot(t,x,'b.');
hold on;
plot(t,x2,'r');
xlim([0 5]);

subplot(223);
plot(lag./Fs,xc);
xlim([-.5 .5]);

subplot(224);
[ax,h1,h2] = plotyy(f,y,f2,(S));
set(ax,'XLim',[0 100]);