t=0:0.001:5;
f=16;
beta=sin(2*pi*t*f);

[loc]=findpeaks(beta,.99);
loc = loc.loc;
spks=zeros(length(t),1);
for n=1:length(loc)
    jitt=floor(rand(1,1)*20)-5;
    jloc(n)=loc(n)+jitt;
end
spks(loc,1)=1;
figure;
subplot(3,1,1);plot(t,beta);
subplot(3,1,2);plot(t,spks);

[A,B]=xcorr(spks);
subplot(3,1,3);plot(B,A);axis([5 250 0 30]);

%%
addpath(genpath('/home/rouxf/tbx/chronux_2_11/spectral_analysis/'));

dt = t(2)-t(1);
Fs = 1/dt;
TW = 3;%(size(data1,1)/Fs)*1/(size(data1,1)/Fs);
params1                  = [];
params1.tapers           = [TW , 2*TW-1];
params1.pad              = -1;
params1.Fs               = 1/dt;
params1.err              = 0;
params1.trialave         = 1;
params1.fpass            = [0 100];

[S,f] = mtspectrumpb(spks,params1);

figure;
plot(f,S);