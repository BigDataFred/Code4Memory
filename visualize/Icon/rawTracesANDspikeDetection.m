
x = [];
for it = 1:length(spkIdx)
    
    if it <2
        x = [x spkIdx{it}];
    else
        x = [x spkIdx{it}+sum([nsamp{it-1:-1:1}])];
    end;
end;
spkIdx2 = x;

Wp = [500 5000]./params.hdr.Fs;
[b,a] = butter(4,Wp);

params.Hd{1} = b;
params.Hd{2} = a;

filteredSignal = filterSignal( params.Hd, dataSamples );

figure;
subplot(211);
plot(dataSamples);
%xlim([1e5 1.5e5]);
subplot(212);
hold on;
plot(filteredSignal,'Color',[.5 .5 .5]);
plot(spkIdx2,filteredSignal(spkIdx2),'r.');
%xlim([1e5 1.5e5]);

%%
[lfp_interp] = interpLFP(dataSamples,spkIdx2,[0.001 0.002],params.hdr2.Fs);