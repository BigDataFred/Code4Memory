%%
selIx7 = find(ttls ==7);
selIx7 = selIx7(trlENC);

dum1 = selIx7(find(ttls(selIx7+1)==0))+1;

dum2 = find(ttls(selIx7+1) ~=0);
dum2 = selIx7(dum2(ttls(selIx7(dum2)+2) ==0))+2;

selIx0 = sort([dum1 dum2]);

%%
selIx = selIx7;

dt = [-0.5*Fs 15*Fs];
%dt = [1*32 15*32];

x = [];
for it =1:length(selIx)    
    x(it,:) = dataSamples(TTLix(selIx(it))-dt(1):TTLix(selIx(it))+dt(2));    
end;

tx = [-dt(1):dt(2)]./Fs;

figure;
plot(tx,mean(x,1),'LineWidth',3);
axis tight;
xlabel('Time rel. to TTL-timestamp [ms]');
ylabel('LFP-Amplitude [\muV]');
box off;

%%
selIx = selIx7;

dt = [32*15 32*15];

x = [];
for it =1:length(selIx)    
    x(it,:) = dataSamples(TTLix(selIx(it))-dt(1):TTLix(selIx(it))+dt(2));    
end;

[b,a] = butter(4,300/(Fs/2),'high');% apply low-pass for LFP
dum = filtfilt(b,a,dataSamples);
x2 = [];
for it = 1:size(x,1)
    [x2(it,:)] = dum(TTLix(selIx(it))-dt(1):TTLix(selIx(it))+dt(2));
end;

tx = [-dt(1):dt(2)]./Fs.*1e3;

figure;
subplot(4,1,1:3);
imagesc(tx,1:size(x,1),x);
caxis([-5 5]);
axis xy;
ylabel('TTL event #');
subplot(4,1,4);
plot(tx,mean(x,1),'LineWidth',3);
axis tight;
xlabel('Time rel. to TTL-timestamp [ms]');
ylabel('LFP-Amplitude [\muV]');

figure;
subplot(4,1,1:3);
imagesc(tx,1:size(x,1),x2);
caxis([-5 5]);
axis xy;
ylabel('TTL event #');
subplot(4,1,4);
plot(tx,mean(x2,1),'LineWidth',3);
axis tight;
xlabel('Time rel. to TTL-timestamp [ms]');
ylabel('LFP-Amplitude [\muV]');

%%
selIx = selIx0;

dt = [32*15 32*15];

x = [];
for it =1:length(selIx)    
    x(it,:) = dataSamples(TTLix(selIx(it))-dt(1):TTLix(selIx(it))+dt(2));    
end;

[b,a] = butter(4,300/(Fs/2),'high');% apply low-pass for LFP
dum = filtfilt(b,a,dataSamples);
x2 = [];
for it = 1:size(x,1)
    [x2(it,:)] = dum(TTLix(selIx(it))-dt(1):TTLix(selIx(it))+dt(2));
end;

tx = [-dt(1):dt(2)]./Fs.*1e3;

figure;
subplot(4,1,1:3);
imagesc(tx,1:size(x,1),x);
caxis([-5 5]);
axis xy;
ylabel('TTL event #');
subplot(4,1,4);
plot(tx,mean(x,1),'LineWidth',3);
axis tight;
xlabel('Time rel. to TTL-timestamp [ms]');
ylabel('LFP-Amplitude [\muV]');

figure;
subplot(4,1,1:3);
imagesc(tx,1:size(x,1),x2);
caxis([-5 5]);
axis xy;
ylabel('TTL event #');
subplot(4,1,4);
plot(tx,mean(x2,1),'LineWidth',3);
axis tight;
xlabel('Time rel. to TTL-timestamp [ms]');
ylabel('LFP-Amplitude [\muV]');

%%
%[b,a] = butter(4,300/Fs,'high');
%dum = filtfilt(b,a,dataSamples);
dum = dataSamples;

selIx = selIx7;

dt = [32*5 32*5];

x = [];
for it =1:length(selIx)    
    x(it,:) = dum(TTLix(selIx(it))-dt(1):TTLix(selIx(it))+dt(2));    
end;

tx = [-dt(1):dt(2)]./Fs.*1e3;

figure;
subplot(4,1,1:3);
imagesc(tx,1:size(x,1),x);
caxis([-5 5]);
axis xy;
ylabel('TTL event #');
subplot(4,1,4);
hold on;
plot(tx,x-repmat(mean(x,2),[1 size(x,2)]),'b');
plot(tx,mean(x,1),'r','LineWidth',3);
axis tight;
xlabel('Time rel. to TTL-timestamp [ms]');
ylabel('LFP-Amplitude [\muV]');
yl = get(gca,'YLim');

dt = [32*0.25 32*2];
[cleanDat,res] = cleanARTIFACTfromLFP(dataSamples,TTLix(selIx),dt,2,2,Fs);
%[b,a] = butter(4,300/Fs,'high');
%cleanDat = filtfilt(b,a,cleanDat);

dt = [32*15 32*15];
tx = [-dt(1):dt(2)]./Fs.*1e3;
x = [];
x2 = [];
for it =1:length(selIx)    
    x(it,:) = cleanDat(TTLix(selIx(it))-dt(1):TTLix(selIx(it))+dt(2));    
end;

figure;
subplot(4,1,1:3);
imagesc(tx,1:size(x,1),x);
caxis([-5 5]);
axis xy;
ylabel('TTL event #');
subplot(4,1,4);
hold on;
plot(tx,x-repmat(mean(x,2),[1 size(x,2)]),'b');
plot(tx,mean(x,1),'r','LineWidth',3);
axis tight;
xlabel('Time rel. to TTL-timestamp [ms]');
ylabel('LFP-Amplitude [\muV]');
ylim(yl);

%%
[b,a] = butter(4,300/Fs,'high');
dum = filtfilt(b,a,dataSamples);
%dum = dataSamples;

selIx = selIx0;

dt = [32*15 32*15];

x = [];
for it =1:length(selIx)    
    x(it,:) = dum(TTLix(selIx(it))-dt(1):TTLix(selIx(it))+dt(2));    
end;

tx = [-dt(1):dt(2)]./Fs.*1e3;

figure;
subplot(4,1,1:3);
imagesc(tx,1:size(x,1),x);
caxis([-5 5]);
axis xy;
ylabel('TTL event #');
subplot(4,1,4);
hold on;
plot(tx,x-repmat(mean(x,2),[1 size(x,2)]),'b');
plot(tx,mean(x,1),'r','LineWidth',3);
axis tight;
xlabel('Time rel. to TTL-timestamp [ms]');
ylabel('LFP-Amplitude [\muV]');
yl = get(gca,'YLim');

dt = [32*1.5 32*3];
[cleanDat] = cleanARTIFACTfromLFP(dataSamples,TTLix(selIx),dt,3,2,Fs);
[b,a] = butter(4,300/Fs,'high');
cleanDat = filtfilt(b,a,cleanDat);

dt = [32*15 32*15];
tx = [-dt(1):dt(2)]./Fs.*1e3;
x = [];
x2 = [];
for it =1:length(selIx)    
    x(it,:) = cleanDat(TTLix(selIx(it))-dt(1):TTLix(selIx(it))+dt(2));    
end;

figure;
subplot(4,1,1:3);
imagesc(tx,1:size(x,1),x);
caxis([-5 5]);
axis xy;
ylabel('TTL event #');
subplot(4,1,4);
hold on;
plot(tx,x-repmat(mean(x,2),[1 size(x,2)]),'b');
plot(tx,mean(x,1),'r','LineWidth',3);
axis tight;
xlabel('Time rel. to TTL-timestamp [ms]');
ylabel('LFP-Amplitude [\muV]');
ylim(yl);

%%
dataSamples2 = dataSamples;%-mean(dataSamples);
%[b,a] = butter(4,300/(Fs/2),'high');% apply low-pass for LFP
%[dataSamples2] = filtfilt(b,a,dataSamples2);

selIx = selIx7;
[cleanDat] = cleanARTIFACTfromLFP(dataSamples2,TTLix(selIx),[32*1 32*3],3,2,Fs);

selIx = selIx0;
[cleanDat] = cleanARTIFACTfromLFP(cleanDat,TTLix(selIx),[32*1 32*3],3,2,Fs);
[cleanDat2] = cleanARTIFACTfromLFP(dataSamples2,TTLix(selIx),[32*1 32*3],3,2,Fs);

% [b,a] = butter(4,300/(Fs/2),'high');% apply low-pass for LFP
% [dataSamples2] = filtfilt(b,a,dataSamples2);
% [cleanDat] = filtfilt(b,a,cleanDat);
% [cleanDat2] = filtfilt(b,a,cleanDat2);

selIx = selIx7;
selIx2 = selIx0;

dt = [1*Fs 5*Fs];

x = [];
x2 = [];
x3 = [];
x4 = [];
for it =1:length(selIx)    
    x(it,:) = dataSamples2(TTLix(selIx(it))-dt(1):TTLix(selIx(it))+dt(2));
    x2(it,:) = cleanDat(TTLix(selIx(it))-dt(1):TTLix(selIx(it))+dt(2));  
    x3(it,:) = dataSamples2(TTLix(selIx2(it))-dt(1):TTLix(selIx2(it))+dt(2));
    x4(it,:) = cleanDat2(TTLix(selIx2(it))-dt(1):TTLix(selIx2(it))+dt(2));
end;

tx = [-dt(1):dt(2)]./Fs;

figure;
subplot(3,2,1:2);
hold on;
plot(tx,mean(x,1),'LineWidth',3);
plot(tx,mean(x2,1),'r');
axis tight;
xlabel('Time rel. to TTL-timestamp [ms]');
ylabel('LFP-Amplitude [\muV]');
subplot(3,2,3);
hold on;
plot(tx,mean(x,1),'LineWidth',3);
plot(tx,mean(x2,1),'r');
axis tight;
xlabel('Time rel. to TTL-timestamp [ms]');
ylabel('LFP-Amplitude [\muV]');
xlim([-.025 .025]);
subplot(3,2,4);
hold on;
plot(tx,mean(x,1),'LineWidth',3);
plot(tx,mean(x2,1),'r');
axis tight;
xlim([1.95 2]);
xlabel('Time rel. to TTL-timestamp [ms]');
ylabel('LFP-Amplitude [\muV]');
subplot(3,2,5);
hold on;
plot(tx,mean(x3,1),'LineWidth',3);
plot(tx,mean(x4,1),'r');
axis tight;
xlim([-.025 .025]);
xlabel('Time rel. to TTL-timestamp [ms]');
ylabel('LFP-Amplitude [\muV]');

%%
movingwin               = [0.25 0.001];
T = movingwin(1);
W = 4;
TW = T*W;
k = 2*TW-1;

params                  = [];
params.pad              = 0;
params.fpass            = [0 30];
params.Fs               = Fs;
params.tapers           = [TW k];
params.trialave         = 0;

[S1,t,f] = mtspecgramc( gradient(x)', movingwin, params );
[S2,t,f] = mtspecgramc( gradient(x2)', movingwin, params );

figure;
subplot(4,1,1:3);
imagesc(t-1,f,10*log10(squeeze(mean(S1,3)))');
axis xy;
xlim([-.4 4]);
subplot(4,1,4);
plot(tx,mean(x,1));
axis tight;
xlim([-.4 4]);

figure;
subplot(4,1,1:3);
imagesc(t-1,f,10*log10(S2)');
axis xy;
xlim([-.4 4]);
subplot(4,1,4);
plot(tx,mean(x2,1));
axis tight;
xlim([-.4 4]);

%%
for jt = 1:length(sortedSpikes)
    [cID] = unique(sortedSpikes{jt}.assignedClusterSeg);
    for it = 1:length(cID)
        
        ix1 = find(sortedSpikes{jt}.assignedCluster == cID(it));
        
        ix = find(sortedSpikes{jt}.assignedClusterSeg == cID(it));
        trl = sortedSpikes{jt}.trl(ix);
        ts = sortedSpikes{jt}.SpikeTimesSeg(ix).*1e3;
        
        trlID = 1:52;%unique(trl);
        if length(ix)>50
            figure;
            subplot(121);
            hold on;
            plot(linspace(0,2,64),sortedSpikes{jt}.wavf(sortedSpikes{jt}.oriIx(ix),:),'r');
            plot(linspace(0,2,64),mean(sortedSpikes{jt}.wavf(sortedSpikes{jt}.oriIx(ix),:),1),'k','LineWidth',3);
            axis tight;
            xlabel('Time (ms)');
            ylabel('Amplitude [\muV]');
            title(['chan:',num2str(jt),' ',chanLab{jt},' cluster:',num2str(cID(it))]);
            
            subplot(122);
            hold on;
            for kt = 1:length( trlID )
                ix2 = find(trl == trlID(kt));
                x = ts(ix2);
                x = [x;x];
                y = kt*ones(1,length(x));
                y = [y-.5;y+.5];
                line(x,y,'Color','k');
                
            end;
            ylim([0 length(trlENC)+1])
            xlim([-500 4e3]);
            xlabel('Time rel. to Cue (ms)');
            ylabel('Trial #');
        end;
    end;
end;
