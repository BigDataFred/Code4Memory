%%
addpath('/home/rouxf/tbx/releaseDec2015/');
addpath(genpath('/home/rouxf/tbx/osort-v3-rel/'));
addpath(genpath('/home/rouxf/tbx/wave_clus/'));
addpath(genpath('/home/rouxf/tbx/chronux_2_11/'));

addpath('/home/rouxf/prj/Bham/code/mcode/params/');
addpath(genpath('/home/rouxf/prj/Bham/code/mcode/helper/'));
addpath('/home/rouxf/prj/Bham/code/mcode/utils/spikes/');

%%
Fs =32e3;
[par] = set_parameters_Bham(Fs);

%%
%[rpath] = '/media/rouxf/rds-share/Archive/MICRO/P09/Spontaneous/';

[rpath] = '/media/rouxf/rds-share/Archive/MICRO/P07/cnEM/';

x = dir(rpath);
x(1:2) = [];

[p2d] = [rpath,x(2).name,'/'];
[CSCfiles] = dir([p2d,'*.ncs']);

[savepath] = p2d;
savepath = [savepath(1:regexp(savepath,'Archive')-1),'iEEG_DATA',savepath(regexp(savepath,'Archive/')+7:end)];
chck = dir(savepath);
if isempty(chck)
    mkdir(savepath);
end;

%%
FieldSelection(1) = 1;%timestamps
FieldSelection(2) = 0;% EventIDs
FieldSelection(3) = 1;%TTLs
FieldSelection(4) = 0;% Extras
FieldSelection(5) = 0;% Event strings

ExtractHeader = 1;

ExtractMode = 1;

ModeArray = [];

[EVfile] = dir([p2d,'*.nev']);

[TimeStamps, ttls, Hdr] = Nlx2MatEV_v3( [p2d,EVfile.name], FieldSelection, ExtractHeader, ExtractMode, ModeArray );

[events] = zeros(size(ttls,2),2);
events(:,1) = TimeStamps';
events(:,2) = ttls';

%% read the CSC data
FieldSelection(1) = 1;%timestamps
FieldSelection(2) = 0;
FieldSelection(3) = 0;%sample freq
FieldSelection(4) = 0;
FieldSelection(5) = 1;%samples
ExtractHeader = 1;

ExtractMode = 1; % 2 = extract record index range; 4 = extract timestamps range.

%%
if isempty(gcp('nocreate'))
    parpool(36,'SpmdEnabled',false);
end;

%%
c = 0;
cluInf = [];
parfor it = 1:length(CSCfiles)
    
    tic;
    fprintf([num2str(it),'/',num2str(length( CSCfiles) )]);
    
    %%
    dum = CSCfiles(it).name;
    dum(regexp(dum,'_')) = [];
    dum(regexp(dum,'CSC'):regexp(dum,'CSC')+2) = [];
    dum(regexp(dum,'.ncs'):end) = [];
    chanLab = dum;
    
    %%
    
    [timestamps, dataSamples,hdr] = Nlx2MatCSC_v3([p2d,CSCfiles(it).name], FieldSelection, ExtractHeader, ExtractMode, []);
    
    chck = regexp(hdr,'ADBitVolts');
    selIdx = [];
    for jt = 1:length(chck);selIdx(jt) = ~isempty(chck{jt});end;
    selIdx = find(selIdx~=0);
    scalef = str2double(hdr{selIdx}(min(regexp(hdr{selIdx},'\d')):end));
    
    
    chck = regexp(hdr,'SamplingFrequency');
    selIdx = [];
    for jt = 1:length(chck);selIdx(jt) = ~isempty(chck{jt});end;
    selIdx = find(selIdx~=0);    
    Fs = str2double(hdr{selIdx}(min(regexp(hdr{selIdx},'\d')):end));
    
    %% flatten
    dataSamples=double(dataSamples(:))';
    dataSamples = dataSamples.*scalef.*1e6;        
    
    %%    
    [b,a] = butter(4,[300]./(Fs/2),'low');% apply low-pass for LFP
    [LFPsig] = filtfilt(b,a,dataSamples);
    
    %%
    [LFPsig,~] = CleanLineNoise( LFPsig ,'Fs', Fs , 'noiseFreq', 50,'windowSize',1);        
    
    %%
    [lfpTime] = 0:1/Fs:(size(LFPsig,2)-1)/Fs;
    
    %%
    [savename] = [chanLab,'_LFPdat.mat'];
    
    readme = {'lfpDat','Fs','lfpTime','chanLab'};
    par_save([[savepath,'lfp_dat/'],savename],[],{LFPsig,Fs,lfpTime,chanLab,readme});

    %%          
    [~,spikeWaveforms,~,spikeTimestamps,~,~] = amp_detect(dataSamples,par);
    
    %%
    waveclus.spikes                     = spikeWaveforms;
    waveclus.index                      = spikeTimestamps;
    
    [dim] = size(spikeWaveforms);
    
    %%
    sortedSpikes.newSpikeTimes     = [];
    sortedSpikes.assignedCluster   = [];
    sortedSpikes.wavf              = [];
    sortedSpikes.num_clus          = [];
    
    if dim(1) > dim(2)% the number of spike events must larger than the number of samples in the waveform
        [sortedSpikes,wltCoeffs] = doSpikeSorting_waveclus( waveclus , par );
    end;
    waveclus = [];
    
    %%
    [savename] = [chanLab,'_SPKdat.mat'];
    
    readme = {'sortedSpikes','wltCoeffs'};
    par_save([[savepath,'spike_dat/'],savename],[],{sortedSpikes,wltCoeffs,readme});
    
    %%        
    fprintf('\n');
    toc;
    
end;

%%
S = cell(1,length(CSCfiles));
S2 = cell(1,length(CSCfiles));
Fs2 = Fs/32;
LFPsig2 = cell(1,length(LFPsig));

for it = 1:length(LFPsig)
    
    fprintf([num2str(it),'/',num2str(length( CSCfiles) )]);
    
    [LFPsig2{it}] = downsample(LFPsig{it},32);   
    
    %%
    movingwin = [1 .01];
    
    params                  = [];
    params.Fs               = Fs2;
    params.pad              = 0;
    params.fpass            = [1 30];
    params.tapers           = [2 3] ;
    params.trialave         = 0;
    
    [S{it},t,f] = mtspecgramc(gradient(LFPsig2{it})', movingwin, params);
    
    params.tapers           = [64 127] ;
    [S2{it},f2] = mtspectrumc(LFPsig2{it}', params);
    
    fprintf('\n');
    
end;

%%
xc = cell(1,length(CSCfiles));
n2 = cell(1,length(CSCfiles));
ts = cell(1,length(CSCfiles));

for it = 1:length(sortedSpikes1)
    
    fprintf([num2str(it),'/',num2str(length( CSCfiles) )]);
    
    [cluID] = unique(sortedSpikes1{it}.assignedCluster);
    
    for ct = 1:length(cluID)
        
        if ~isempty(cluID) && cluID(ct) ~=0            
            if ~isempty(sortedSpikes1{it}.newSpikeTimes)
                
                c = c+1;
                
                cluInf(c,1) = it;
                cluInf(c,2) = cluID(ct);
                
                ix = find( sortedSpikes1{it}.assignedCluster == cluID(ct) );
                
                ts{c} = lfpTime( sortedSpikes1{it}.newSpikeTimes( ix ) );
                ts{c} = ts{c}.*1e3;
                
                wvf{c} = sortedSpikes1{it}.wavf(ix,:);
                
                isi = diff(ts{c});
                
                [n,~] = hist(ts{c},[0:1:max(ts{c})]);
                [n2{c},bx2] = hist(isi,[-2:1:251]);
                n2{c}(end) = [];bx2(end) = [];
                
                [xc{c},lag] = xcorr(n,250);
                
            end;            
        end;
        
    end;
    
    fprintf('\n');
    
end;

%%
%parpool(18,'SpmdEnabled',false);

lfpTime2 = 0:1/Fs2:(size(LFPsig2{1},2)-1)/Fs2;

[dt] = lfpTime2(1)*1e3:1:lfpTime2(end)*1e3;

S3 = cell(1,length(ts));
R = cell(1,length(ts));
C = cell(length(ts),1);

for it = 1:size(cluInf,1)
    
    fprintf([num2str(it),'/',num2str(length(ts))]);
    
    n = hist(ts{it}',dt);
    
    params                  = [];
    params.Fs               = 1e3;
    params.pad              = 0;
    params.fpass            = [1 30];    
    params.trialave         = 0;        
    params.tapers           = [64 127] ;
    
    [S3{it},f3,R{it}] = mtspectrumpb(n, params);
    
    c = zeros(length(LFPsig2),length(f3));
   for jt = 1:length(LFPsig2)
        [c(jt,:)] = coherencycpb(LFPsig2{jt}',n', params);
    end;
    C{it} = c;
    
    fprintf('\n');
    
end;

%%

params                  = [];
params.Fs               = 1e3;
params.pad              = 0;
params.fpass            = [1 30];
params.trialave         = 0;
params.tapers           = [64 127] ;

[C2,~,~,~,~,f3] = coherencyc(LFPsig2{32}',LFPsig2{43}',params);
[C3,~,~,~,~,f3] = coherencyc(LFPsig2{32}',LFPsig2{45}',params);
[C4,~,~,~,~,f3] = coherencyc(LFPsig2{43}',LFPsig2{45}',params);

%%

dt = lfpTime2.*1e3;
[n1,~] = hist(ts{43},dt);
[n2,~] = hist(ts{45},dt);

movingwin = [1 .01];

params                  = [];
params.Fs               = Fs2;
params.pad              = 0;
params.fpass            = [1 30];
params.tapers           = [2 3] ;
params.trialave         = 0;

[C5,~,~,~,~,t,f4] = cohgramcpb(LFPsig2{32}', n1',movingwin, params);
[C6,~,~,~,~,t,f4] = cohgramcpb(LFPsig2{32}', n2',movingwin, params);
[C7,~,~,~,~,t,f4] = cohgrampb(n1', n2',movingwin, params);

%%
% y= [];
% for it = 37%1:length(C)
%     
%     y(it) = mean(C{it}(f3 >=4 & f3 <=6));
% end;
% [~,ix] = max(y);
ix = 38;

figure;
subplot(131);
y = S2{ix};
x = f2';
y = log10(y);
x = log10(x);
[b] = regress(y,[ones(size(x)) x x.^2 x.^3]);
yp = [ones(size(x)) x x.^2 x.^3]*b;
hold on;
plot(f2,y);
plot(f2,yp,'r');
axis tight;axis square;
subplot(132);
plot(f3,S3{45});
axis tight;axis square;
subplot(133);
plot(f3,C{45});
axis tight;axis square;

%%
dt = lfpTime2(1):1:lfpTime2(end).*1e3;
[n,~] = hist(ts{43},dt);
ix = find(n==1);

c = 0;
STA = [];
for it = 1:length(ix)
    
    if ( sign(ix(it)-1000)==1 ) && ( ix(it)+1000 < length(LFPsig2{1}))
        
        x = LFPsig2{37}(ix(it)-1000:ix(it)+1000);
        c = c+1;
        STA(c,:) = x;
    end;
    
end;

tmp = mean(STA,1);

phi1 = angle(hilbert(tmp));

phi2 = angle(hilbert((STA)')');

plv = 1/size(phi1,2).*abs(sum(exp(1i.*(ones(size(phi2,1),1)*phi1-phi2)),2));
[~,ix] = sort(plv);

figure;
subplot(4,1,1:3);
imagesc(-1000:1000,1:c,(STA(ix,:)));
axis xy;
subplot(4,1,4);
hold on;
%plot(-500:500,STA,'k');
plot(-1000:1000,mean(STA,1),'r');
axis tight;

params                  = [];
params.Fs               = 1e3;
params.pad              = 0;
params.fpass            = [1 30];
params.trialave         = 1;
params.tapers           = [2 3] ;

[S4,f4] = mtspectrumc(mean(STA,1)', params);
[S5,f5] = mtspectrumc(gradient(STA)', params);

movingwin = [.5 .01];
params                  = [];
params.Fs               = 1e3;
params.pad              = 0;
params.fpass            = [1 30];
params.trialave         = 1;
params.tapers           = [2 3] ;

[S6,t6,f6] = mtspecgramc((STA-(ones(size(STA,1),1)*mean(STA,1)))',movingwin,params);

figure;
subplot(121);
hold on;
plot(f4,log10(S4),'b');
plot(f5,log10(S5),'r');
axis square;
subplot(122);
plot(f4,S4./S5);
axis square;

%%
m = [];
a = [];
figure;
for it = 1:length(S3)
    subplot(9,8,it);
    a(it) = gca;
    %plot(f3,S3{it}./R{it});
    plot(S3{it});
    axis tight;
    %m(it,:)= [min(S3{it}./R{it}) max(S3{it}./R{it})];
    m(it,:)= [min(S3{it}) max(S3{it})];
end;
set(a,'Ylim',[min(min(m)) max(max(m))]);

%%
lfpTime2 = 0:1/Fs2:(size(LFPsig2{32},2)-1)/Fs2;
c = [];
for it = 1:length( LFPsig2 )
    
    %figure;      
    y = S2{it};
    x = f2';
    y = log10(y);
    x = log10(x);
    [b] = regress(y,[ones(size(x)) x x.^2 x.^3]);
    yp = [ones(size(x)) x x.^2 x.^3]*b;
    c(it,:) = (y-yp);
       
%     subplot(3,10,1:2);
%     hold on;
%     plot(f2,log10(S2{it}));
%     plot(f2,yp,'r');
%     xlabel('Frequency (Hz)');
%     ylabel('Power (log)');
%     legend('LFP','1/f');
%     axis square;
%     
%     axis tight;
%     subplot(3,10,4:5);
%     plot(f2,c(it,:));
%     axis tight;
%     xlabel('Frequency (Hz)');
%     ylabel('1/f-residual (log)');
%     axis square;
    
%     init = t(1);
%     while init < t(end)
%         subplot(3,10,21:30);
%         ix = 1:length(t);%find(t>=init & t<=init+2);
%         imagesc(t(ix),f,S{it}(ix,:)');
%         axis xy;caxis([0 1]);
%         colormap jet;
%         xlabel('Time (s)');
%         ylabel('Frequency (Hz)');
%         %set(gca,'XTick',round([min(t(ix)):.2:max(t(ix))].*10)/10);
%         
%         subplot(3,10,11:20);
%         ix2 = find(lfpTime2 >= t(ix(1)) & lfpTime2 <= t(ix(end)));
%         %hold on;
%         %plot([min(lfpTime2(ix2)) max(lfpTime2(ix2))],[0 0],'r');
%         plot(lfpTime2(ix2),(LFPsig2{it}(ix2)),'k');
%         axis tight;
%         xlabel('Time (s)');
%         ylabel('Amplitude (\muV)');
%         ylim([-500 500]);
%         %set(gca,'XTick',round([min(lfpTime2(ix2)):.2:max(lfpTime2(ix2))].*10)./10);
%         
%         %pause;
%         init = init+2;
%         
%     end;
    
end;

y= [];
for it = 1:size(c,1)
    
    y(it) = max(c(it,:));    

end;

%%
hemLab = cell(1,length(chanLab));
ix = regexp(chanLab,'\d{1}');
for it = 1:length(ix)    
   ix2 = ix{it}-1;   
   hemLab(it) = {chanLab{it}(ix2)};
end;

dum = {};
for it = 1:size(cluInf,1)
    dum(it) = hemLab(cluInf(it,1));
end;

ix = [];
ix(:,1) = [find(strcmp(dum,'L'))';find(strcmp(dum,'R'))'];
ix(:,2) = [ones(length(find(strcmp(dum,'L'))),1);2*ones(length(find(strcmp(dum,'R'))),1)];

figure;
subplot(5,1,1);
hold on;
% x1 = zeros(1,length(LFPsig2{1}));
% x2 = zeros(1,length(LFPsig2{1}));
% c1 = 0;
% c2 = 0;
% for it  =1:length(LFPsig2)
%     if strcmp(hemLab(it),'L')
%         c1 = c1+1;
%         x1 = x1+LFPsig2{it};
%     else
%         c2 = c2+1;
%         x2 = x2+LFPsig2{it};
%     end;
% end;
% x1 = x1./c1;
% x2 = x2./c2;
%plot(lfpTime2,x1,'r');
%plot(lfpTime2,x2,'b');
%plot(lfpTime2,LFPsig2{ix1},'r');
sel = find(ix(:,2)==2);
for it = 1:length(sel)
plot(lfpTime2,LFPsig2{cluInf(ix(sel(it),1),1)},'k');
end;
axis tight;

subplot(5,1,2:5);
hold on;
for it = 1:length(sel)
    
    x = ts{ix(sel(it),1)}./1e3;
    x = [x;x];
    
    y = it*ones(1,length(x));
    y = [y-.5;y+.5];
    
    if ix(sel(it),2) == 1
        line(x,y,'Color','r');
    else
        line(x,y,'Color','b');
    end;
    
end;

%%
for it = 43;%ix%1:size(cluInf,1)
    figure;
    subplot(2,2,1);
    hold on;
    plot(linspace(0,2,64),wvf{it},'r');
    plot(linspace(0,2,64),mean(wvf{it},1),'k','LineWidth',3);
    plot([1.5 2],[min(min(wvf{it})) min(min(wvf{it}))],'k','LineWidth',3);
    plot([2 2],[min(min(wvf{it})) min(min(wvf{it}))+75],'k','LineWidth',3);
    axis off;axis square;
    title([chanLab{cluInf(it,1)},' ,cluster #',num2str(cluInf(1,2))]);
    subplot(2,2,3);
    bar(bx2,n2{it});axis tight;axis square;
    subplot(2,2,4);
    selIx = find(lag==0)+1:length(lag);
    bar(lag(selIx),xc{it}(selIx));
    axis tight;axis square;
end;

%%
dt = min(lfpTime2):1:max(lfpTime2).*1e3;
[n1,~] = hist(ts{43},dt);
[n2,~] = hist(ts{45},dt);

[xc,lag] = xcorr(n1,n2,500);


[b,a] = butter(4,[10]./(Fs2/2),'low');% apply low-pass for LFP
filtXC = filtfilt(b,a,xc);

figure;
hold on;
bar(lag,xc);
plot([0 0],[0 max(xc)],'Color',[.75 .75 .75],'LineWidth',3);
plot(lag,filtXC,'r','LineWidth',3);
axis tight;
set(gca,'XTick',[-500 -250 0 250 500]);
axis tight; axis square;


%%
lag = 500;
xc = zeros(lag*2+1,length(dt));

parfor kt = 1:length( dt ) % loop over bins
    
    dum1 = n1;dum2 = n2;
    
    if (kt>lag) && (kt < length( dum1 )-lag)
        ix2 = kt-lag/2:kt+lag/2;% shift window for autocorrelation in 1 ms steps
        dum = xcorr(dum1(ix2),dum2(ix2),lag);% compute autocorrelation
        xc(:,kt) = dum';% sum autocorr across trials for each bin
    end;
    
end;

%%
figure;
hold on;
x =  ts{43};
y = ones(1,length(x));
x = [x;x];
y = [y-.5;y+.5];
line(x./1e3,y,'Color','k');

x =  ts{45};
y = 2*ones(1,length(x));
x = [x;x];
y = [y-.5;y+.5];
line(x./1e3,y,'Color','k');

%%
[chanG] = {};
for it = 1:length(chanLab)
    
    chanG(it) = {chanLab{it}(1:regexp(chanLab{it},'\d{1}')-1)};
    
end;

[chanID] = unique( chanG );

avgLFP = {};
for it = 1:length( chanID )
    [ix] = find(strcmp(chanG,chanID(it)));
    
    n = length(ix);
    avgLFP{it} = zeros(1,length(lfpTime2));
    for jt = 1:length( ix )
        avgLFP{it} = avgLFP{it}+LFPsig2{ix(jt)};
    end;
    avgLFP{it} = avgLFP{it}./n;
    
end;

figure;
for it = 1:length(avgLFP)
    subplot(length(avgLFP),1,it);
    plot(lfpTime2,avgLFP{it});
    axis tight;
end;

params                  = [];
params.Fs               = Fs2;
params.pad              = 0;
params.fpass            = [1 30];
params.trialave         = 0;
params.tapers           = [56 127] ;

n = length(avgLFP);
np = 0;
for it = 1:n-1
    np = np+n-it;    
end;

lfpCoh = cell(np,1);
ix = 1:length(avgLFP);
k = 0;
for it = 1:length(avgLFP)
    fprintf([num2str(it),'/',num2str(length(avgLFP))]);
    ix2 = setdiff(ix,ix(1));
    for jt = 1:length(ix2)
        k = k+1;
        [lfpCoh{k},~,~,~,~,f3] = coherencyc(avgLFP{ix(1)}',avgLFP{ix2(jt)}',params);
    end;
    ix(1) = [];
    fprintf('\n');
end;

k=0;
c = 0;
figure;
ix = 1:n;
for it = 1:n
    ix2 = setdiff(ix,it);
    for jt = 1:length(ix2)
        k = k+1;
        c = c+1;
        subplot(n-1,n-1,c);
        plot(f3,lfpCoh{k});
        axis tight;
        ylim([0 1]);
        title(chanID(ix2(jt)));
        if jt ==1
            ylabel(chanID(ix(1)));
        end;
    end;
    c=c+ix(1);
    ix(1) = [];
end;

%%
params                  = [];
params.Fs               = Fs2;
params.pad              = 0;
params.fpass            = [1 30];
params.trialave         = 0;
params.tapers           = [56 127] ;

[Slfp1,~] = mtspectrumc( avgLFP{7}', params );
[Slfp2,~] = mtspectrumc( avgLFP{8}', params );
[Slfp3,f3] = mtspectrumc( avgLFP{4}', params );

%%
movingwin = [5 .05];

params                  = [];
params.Fs               = Fs2;
params.pad              = 0;
params.fpass            = [1 30];
params.tapers           = [10 19] ;
params.trialave         = 0;

[lfpCoh1,~,~,~,~,f,f4] = cohgramc((avgLFP{7}'-mean(avgLFP{7})),(avgLFP{8}'-mean(avgLFP{8})),movingwin, params );
[lfpCoh2,~,~,~,~,f,f4] = cohgramc((avgLFP{7}'-mean(avgLFP{7})),(avgLFP{4}'-mean(avgLFP{4})),movingwin, params );

%% 46s is money
idx = 45*32e3:46*32e3;%1:32e4;%46+37*32e3:46+38*32e3;%
figure;
set(gcf,'Color','w');cnt = 0;
while idx(end) < length(dataSamples)
    cnt = cnt+1;
    hold on;
    plot(1:length(idx),x(idx)-mean(x(idx)),'b','LineWidth',1.5);
    xsc = length(idx)-(length(idx)/10):length(idx);
    plot(xsc,1.5*min(x(idx))*ones(1,length(xsc)),'k','LineWidth',3);
    axis tight;
    axis off;title(cnt);
    idx = idx+32e3;
    pause;
    clf;
end;

%% 46s is money
idx = 46+37*32e3:46+38*32e3;%45*32e3:46*32e3;
figure;
set(gcf,'Color','w');cnt = 0;
subplot(3,1,1);
hold on;
plot(1:length(idx),x(idx)-mean(x(idx)),'b','LineWidth',1.5);
xsc = length(idx)-(length(idx)/10):length(idx);
plot(xsc,min(x(idx))*ones(1,length(xsc)),'k','LineWidth',3);
plot(xsc(end)*ones(1,2),[min(x(idx)) min(x(idx))+75],'k','LineWidth',3);
axis tight;
axis off;

subplot(3,1,2);
hold on;
plot(1:length(idx),x2(idx)-mean(x2(idx)),'Color',[.75 .75 .75],'LineWidth',1.5);
z = (x2(idx)-mean(x2(idx)))./std(x2(idx));
thr = z > 4;
ix = find(thr ==1);
ix(diff(ix)==1) = [];
for it = 1:length(ix)
    plot(ix(it),min(x2(idx)-mean(x2(idx)))-5,'r^','MarkerFaceColor','r');
end;
axis tight;
axis off;

subplot(3,1,3);
plot(1:length(idx),x3(idx)-mean(x3(idx)),'k','LineWidth',3);
axis tight;
axis off;

%%
selIdx = find(ismember(spikeTimestamps1,idx));
wvf = spikeWaveforms1(selIdx,:);
figure;
hold on;
plot(linspace(0,2,64),wvf,'r');
plot([1.5 2],[min(min(wvf)) min(min(wvf))],'k','LineWidth',3);
plot([2 2],[min(min(wvf)) min(min(wvf))+75],'k','LineWidth',3);
axis off;

%%
set(gca,'Color','none');
savepath = '/home/rouxf/tmp/figuresIcon/';
set(gcf,'PaperPositionMode','auto');
print(gcf,'-r800','-dtiff',[savepath,'spkikeDetection2.tif']);