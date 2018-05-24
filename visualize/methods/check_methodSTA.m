%%
Fs = 1e3;
t = 0:1/Fs:10-1/Fs;
lfp = sin(2*pi*5.*t);

ix = find(lfp == max(lfp));

spk = zeros(1,length(lfp));
spk(ix) = 1;

ts = t(ix);

dt = min(t):1/Fs:max(t);
raster = hist(ts,dt);

figure;
subplot(311);
plot(t(1:1e3),lfp(1:1e3));axis tight;
subplot(312);
plot(t(1:1e3),spk(1:1e3),'k');axis tight;
subplot(313);
plot(dt(1:1e3),raster(1:1e3),'bo');axis tight;

sig = [];
ix = find(raster == 1);
dt = 250;c=0;
for it = 1:length(ix)
    start = ix(it)-dt;
    stop = ix(it)+dt;
    if (start >0) && (stop < length(lfp))
        c = c+1;
        sig(c,:) = lfp(start:stop);
    end;
end;

figure;
subplot(4,1,1:3);
imagesc(-dt:dt,1:size(sig,1),sig);
subplot(4,1,4);
plot(-dt:dt,mean(sig,1));
axis tight;
ylim([-1 2]);