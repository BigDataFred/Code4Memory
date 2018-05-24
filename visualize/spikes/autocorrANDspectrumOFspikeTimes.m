%%
addpath(genpath('/media/rouxf/rds-share/Common/chronux_2_11/'));

%%
trl = spike_data.trial{tinf{sel(jt)}(1)};
trlID = unique(trl);

ts = spike_data.timestamp{tinf{sel(jt)}(1)};
ts = ts./1e6;
ts = ts.*1e3;
ts = ts-min(ts);
ts = ts -mean(ts);
dt = min(ts):1:max(ts);
n = histc(ts,dt);

[xc,lag] = xcorr(n,100,'coeff');

x = spike_data.time{tinf{sel(jt)}(1)}.*1e3;

r = {};
for it = 1:length(trlID)
    ix = find(trl == trlID(it));
    r{it} = x(ix);
end;

%%
dt = min(x):max(x);
n = zeros(length(dt),length(r));
for it = 1:length(r)
    n(:,it) = histc(r{it},dt);
end;

%%
params = [];
params.Fs = 1e3;
params.fpass = [0 20];
params.pad =1;
params.tapers = [1.5 2];
params.err = [1 0.05];
params.trialave = 0;

[S1,f1,R1,Serr1] = mtspectrumpb(n,params,0);

params = [];
params.Fs = 1e3;
params.fpass = [20 100];
params.pad =1;
params.tapers = [7.5 14];
params.err = [1 0.05];
params.trialave = 0;

[S2,f2,R2,Serr2] = mtspectrumpb(n,params,0);

%%
R1 = R1*ones(1,length(f1));
R1 = R1';
R2 = R2*ones(1,length(f2));
R2 = R2';

%%

figure;
subplot(121);
hold on;
plot(f1,mean(S1./R1,2));
plot(f1,mean(squeeze(Serr1(1,:,:))./R1,2),'r');
plot(f1,mean(squeeze(Serr1(2,:,:))./R1,2),'r');
xlabel('Frequency (Hz)');

subplot(122);
hold on;
plot(f2,mean(S2./R2,2));
plot(f2,mean(squeeze(Serr2(1,:,:))./R2,2),'r');
plot(f2,mean(squeeze(Serr2(2,:,:))./R2,2),'r');
xlabel('Frequency (Hz)');