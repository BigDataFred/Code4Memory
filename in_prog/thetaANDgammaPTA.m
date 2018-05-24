
%% LFP-power and SFC specrta
tInt = [0 1;1 2;2 3;3 4;4 5];

for kt = 1:length(data_lfp)
    
    if any(foi(pvalPL{kt}<0.05/6) >=3) && any(foi(pvalPL{kt}<0.05/6) <=8)
        
        selIdx1 = find(MWdat.time{1} >= -1 & MWdat.time{1} <0); % extract time CUe
        [spectrumLFP1,spectrumSFC1] = computeSpectrumAndSFC4timeRange(data_lfp(kt), data_spk(kt), selIdx1,Fs);
        
        AVG1 = zeros(size(spectrumLFP1{1}.S1,1),1);
        AVG2 = zeros(size(spectrumLFP1{1}.S2,1),1);
        AVG3 = zeros(size(spectrumSFC1{1}.SFC1,1),1);
        AVG4 = zeros(size(spectrumSFC1{1}.SFC2,1),1);
        for it = 1:length(tInt)
            selIdx2 = find(MWdat.time{1} >=tInt(it,1) & MWdat.time{1} <tInt(it,2));% extract time Encoding
            [Pyy,SFCyy] = computeSpectrumAndSFC4timeRange(data_lfp(kt), data_spk(kt), selIdx2,Fs);
            AVG1 =AVG1+mean(Pyy{1}.S1,2);
            AVG2 =AVG2+mean(Pyy{1}.S2,2);
            AVG3 =AVG3+mean(SFCyy{1}.SFC1,2);
            AVG4 =AVG4+mean(SFCyy{1}.SFC2,2);
        end;
        spectrumLFP2 = Pyy;
        spectrumSFC2 = SFCyy;
        
        spectrumLFP2{1}.S1 = AVG1./it;
        spectrumLFP2{1}.S2 = AVG2./it;
        spectrumSFC2{1}.SFC1 = AVG3./it;
        spectrumSFC2{1}.SFC2 = AVG4./it;
        
        [cleanSpectrum1,freqBinz1] = oneOverfRegression(spectrumLFP1{1}.fx1,mean(spectrumLFP1{1}.S1,2));
        [cleanSpectrum2,freqBinz2] = oneOverfRegression(spectrumLFP2{1}.fx1,mean(spectrumLFP2{1}.S1,2));
        
        %%
        smp =MWdat.time{1};
        D = [-.5 .5];
        plt = 'n';
        err = 1;
        
        T = [0 2];
        [s1,t,E] = sta(data_spk2{kt},data_lfp{kt},smp,plt,[],T,D,err);
        
        T = [2.5 4.5];
        [s2,t,E] = sta(data_spk2{kt},data_lfp{kt},smp,plt,[],T,D,err);
        
        dum1 = {mean(s1,1)'};
        dum2 = {mean(s2,1)'};
        
        [spectrumSTA1] = computeSpectrumAndSFC4timeRange(dum1, [], [],Fs);
        [spectrumSTA2] = computeSpectrumAndSFC4timeRange(dum2, [], [],Fs);
        
        [cleanSpectrum3,freqBinz3] = oneOverfRegression(spectrumSTA1{1}.fx1,mean(spectrumSTA1{1}.S1,2));
        [cleanSpectrum4,freqBinz4] = oneOverfRegression(spectrumSTA2{1}.fx1,mean(spectrumSTA2{1}.S1,2));
        
        % %demean and detrend
        % cfg                     = [];
        % cfg.bpfilter            = 'yes';
        % cfg.bpfreq              = [2.75 60];
        % cfg.bpfilttype          = 'but';
        % cfg.bpfiltord           = 4;
        %
        % [dum] = preprocLFP(MWdat , cfg );
        % [data_lfpB] = ft2chronux([],dum, BFinf, trlENC,[min(MWdat.time{1}) max(MWdat.time{1})]);
        %
        % T = [0 2];
        % [s1,t,E1] = sta(data_spk2{1},data_lfpB{1},smp,plt,[],T,D,err);
        %
        % T = [2.5 4.5];
        % [s2,t,E2] = sta(data_spk2{1},data_lfpB{1},smp,plt,[],T,D,err);
        
        visualizeSFCandSTAandSTP;
        
    end;
end;

%%
selIdx1 = find(MWdat.time{1} >= 0 & MWdat.time{1} <2); % extract time CUe
selIdx2 = find(MWdat.time{1} >=2.75 & MWdat.time{1} <4.75);% extract time Encoding

ix1 = [];
for it = 4:6
    ix1 = [ix1;logDat.ix{it}];
end;

ix2 = [];
for it = 4:6
    ix2 = [ix2;logDat.ix{it}];
end;

dz = cell(1,length(data_spk));
vdz = cell(1,length(data_spk));
Adz = cell(1,length(data_spk));
cnt = 0;
for it = 1
    cnt = cnt+1;
    data1c1 = data_lfp{5}(selIdx2,ix1);% s2';%
    data1c2 = data_lfp{5}(selIdx1,ix2);% s1';%
    
    data2c1 = data_spk{3}(selIdx2,ix1);
    data2c2 = data_spk{3}(selIdx1,ix2);
    
    TW = 3;
    
    params                  = [];
    params.Fs               = Fs;
    params.pad              = 2;
    params.fpass            = [0 30];
    params.tapers           = [TW 2*TW-1];
    
    [dz{cnt},vdz{cnt},Adz{cnt},fBinz] = runAversenUtest(data1c1,data1c2,data2c1,data2c2,params,5e-2,1);
    if any(Adz{cnt} ~=1);
        figure;plot(fBinz,dz{cnt});
    end;
end;
   
%%
%[trXC2.xc,trXC2.fr,trXC2.lag] = time_resolved_xcorr([-1e3:5e3],250,trlENC,spkDat{4},spkDat{2});

%%
[cleanSpectrum1,freqBinz] = oneOverfRegression(spectrumSTA{1}.fx1,mean(spectrumSTA{1}.S1,2));
[cleanSpectrum2,freqBinz] = oneOverfRegression(spectrumLFP4{1}.fx1,mean(spectrumLFP4{1}.S1,2));

figure;
hold on;
plot(freqBinz,cleanSpectrum1);
plot(freqBinz,cleanSpectrum2,'r');

%%
figure;
for it = 1:length(spectrumLFP1)
    subplot(1,length(spectrumLFP1),it);
    hold on;
    %plot(spectrumLFP1{it}.fx1,mean(spectrumLFP1{it}.S1,2));
    plot(spectrumLFP1{it}.fx1,spectrumLFP4{it}.S1,'r');
    axis tight;
end;
figure;
for it = 10%:length(spectrumSFC1)
    %subplot(1,length(spectrumSFC1),it);
    %hold on;
    figure;
    hold on;
    %plot(spectrumSFC1{1,it}.fx1,spectrumSFC1{1,it}.SFC1);
    plot(spectrumSFC2{1,it}.fx1,spectrumSFC2{1,it}.SFC1','r');
    axis tight;
end;
figure;
for it = 1:length(spectrumLFP1)
    hold on;
    plot(spectrumLFP1{it}.fx2,mean(spectrumLFP1{it}.S2,2));
    plot(spectrumLFP1{it}.fx2,mean(spectrumLFP2{it}.S2,2),'r');
    axis tight;
end;

for it = 1:length(spectralData.medSgmLFP)
    figure;
    imagesc(spectralData.tx-2,spectralData.fx,10*log10(spectralData.medSgmLFP{it})');
    axis xy;
end;

for it = 1:length(spectralData.Cgrm)
    figure;
    imagesc(spectralData.tx-2,spectralData.fx,spectralData.Cgrm{it}');
    axis xy;caxis([0 .2]);
end;

%%
n = 2; 
n2 = length(spkDat);


figure;
cnt = 0;
cnt2 = length(SgmSPK);
for it = 1:length(SgmSPK)
    cnt = cnt+1;
    subplot(n,n2,cnt);
    imagesc(tx-2,fx2,(SgmSPK2{it}./R2{it}(:,1:length(fx2)))');
    axis xy;
    cnt2 = cnt2+1;
    subplot(n,n2,cnt2);
    imagesc(tx-2,fx,(SgmSPK{it}(:,1:length(fx)))');
    axis xy;
end;

n = length(chID)/8; 
if n >1
    n2 = 8;
else
    n2 = 1;
end;

figure;
for it = 1:length(medSgmLFP)    
    subplot(n,n2,it);
    imagesc(tx-2,fx,10*log10(medSgmLFP{it})');
    axis xy;            
end;

for it = 1:size(Cgrm,1)
    figure;
    cnt = 0;
    cnt2 =size(Cgrm,2);
    for jt = 1:size(Cgrm,2)
                
        cnt = cnt+1;
        subplot(2,size(Cgrm,2),cnt);
        imagesc(tx-2,fx,Cgrm2{it,jt}');
        axis xy;
        
        cnt2 = cnt2+1;
        subplot(2,size(Cgrm,2),cnt2);
        imagesc(tx-2,fx,Cgrm{it,jt}');
        axis xy;
        
    end;
end;

figure;
cnt = 0;
for it = 1:size(Cgrm3,1)
    for jt = 1:size(Cgrm3,2)
        cnt = cnt+1;
        subplot(size(Cgrm3,1),size(Cgrm3,2),cnt);
        imagesc(tx-2,fx,Cgrm3{it,jt}');
        axis xy;
    end;
end;

figure;
cnt = 0;
for it = 1:size(Cgrm4,1)
    for jt = 1:size(Cgrm4,2)
        cnt = cnt+1;
        subplot(size(Cgrm4,1),size(Cgrm4,2),cnt);
        imagesc(tx-2,fx,Cgrm4{it,jt}');
        axis xy;
    end;
end;

%%
[dt] = MWdat.time{1}(2)-MWdat.time{1}(1);
Fs = 1/dt;

TW = 3;
params                  = [];
params.tapers           = [TW , 2*TW-1];
params.pad              = -1;
params.Fs               = Fs;
params.err              = 0;
params.trialave         = 1;
params.fpass            = [0 30];

if isempty(gcp('NoCreate'))
    parpool(36,'SpmdEnabled',false);
end;

n1 = length(data_lfp);
n2 = length(data_spk);
Cemp  = cell(n1,n2);
Cshuf = cell(n1,n2);
pval  = cell(n1,n2);
cnt = 0;
for it = 1:length(medLFPsig)    
    for jt = 1:length(data_spk)        
        cnt = cnt+1;
        [pval{cnt},Cemp{cnt},Cshuf{cnt}] = shuffledSFC(medLFPsig{it},data_spk{jt},params,Fs);
    end;
end;

%%
cnt = 0;
figure;
for it = 1:size(CgmLFPSPK,1)
    for jt = 1:size(CgmLFPSPK,2)
        cnt = cnt+1;
        subplot(8,4,cnt);
        imagesc(tx,fx2,CgmLFPSPK{it,jt}');
        axis xy;
        
    end;
end;

%%
movingwin1 = [0.25 .01];
TW = 2.5;%
params3                  = [];
params3.tapers           = [TW , 2*TW-1];
params3.pad              = -1;
params3.Fs               = Fs;
params3.err              = 0;
params3.trialave         = 0;
params3.fpass            = [30 200];

movingwin2 = [0.25 .01];
TW = 2.5;%
params3                  = [];
params3.tapers           = [TW , 2*TW-1];
params3.pad              = -1;
params3.Fs               = Fs;
params3.err              = 0;
params3.trialave         = 0;
params3.fpass            = [30 200];

SgmLFP = {};
SgmLFP2 = {};

        [   SgmSPK   ] = mtspecgrampb(data_spk{jt},movingwin1,params1);
        
        figure;
        imagesc(tx-1,fx,10*log10(squeeze(mean(SgmSPK,3)))');axis xy;
        title([chID(it) num2str(jt)]);

        figure;
        imagesc(tx-1,fx,CgmLFPSPK');axis xy;
        title([chID(it) num2str(jt)]);

    [SgmLFP{it},tx,fx] = mtspecgramc(data_lfp{1},movingwin1,params1);
    [SgmLFP2{it},~,fx2] = mtspecgramc(data_lfp{1},movingwin3,params3);

figure;
imagesc(tx-1,fx,10*log10(squeeze(mean(SgmLFP{1},3)))');axis xy;

title(chID(it));
figure;
imagesc(tx-1,fx2,10*log10(squeeze(mean(SgmLFP2{1},3)))');axis xy;
title(chID(it));

%%
cfg                     = [];
cfg.preproc             = 'yes';
cfg.preproc             = 'detrend';
cfg.method              = 'wavelet';
cfg.output              = 'fourier';
cfg.pad                 = 'nextpow2';
cfg.taper               ='dpss';
cfg.foi                 = 0:30;
cfg.tapsmofrq           = 1;
cfg.toi                 = dum.time{1}(1):0.01:dum.time{1}(end);
cfg.t_ftimwin           = ones(1,length(cfg.foi));
cfg.keeptrials          = 'yes';

[phi] = ft_freqanalysis(cfg,dum);

itc = phi.fourierspctrm./abs(phi.fourierspctrm);
itc = abs(sum(itc,1)./size(itc,1));
itc = squeeze(itc);

%% bandpass filter the median LFP
cfg                     = [];
cfg.channel             = spkDat{1}.hdr.label;

dum = ft_selectdata( cfg,   MWdat );

cfg                     = [];
%cfg.channel             = spkDat{2}.hdr.label;
cfg.lpfilter            ='yes';
cfg.lpfreq              = 10;
cfg.lpfilttype          = 'but';
%cfg.bpfiltord           = 2;
cfg.padding             = 120;
cfg.padtype             = 'data';

[nbfilt] = ft_preprocessing( cfg , dum );

cfg                     = [];
%cfg.channel             = spkDat{2}.hdr.label;
cfg.bpfilter            ='yes';
cfg.bpfreq              = [1 60];
cfg.bpfilttype          = 'but';
cfg.padding             = 120;
cfg.padtype             = 'data';

[bbfilt] = ft_preprocessing( cfg , dum );

%%
x = nbfilt.trial{1};
ixp = find(sign(diff(sign(diff(x))))==-1);
ixn = find(sign(diff(sign(diff(x))))== 1);
figure;
hold on;
plot(x);
plot(ixp,x(ixp),'r.');
plot(ixn,x(ixn),'g.');

%%
tx2 = -1:1/Fs:4.999;
mix1 = cell(1,length(bbfilt.trial));
mix2 = cell(1,length(bbfilt.trial));

for it = 1:length(bbfilt.trial)
    
    x = bbfilt.trial{it};
    ixp = find(sign(diff(sign(diff(x))))==-1);
    ixn = find(sign(diff(sign(diff(x))))== 1);
    
    tx2(ixp(1))

end;

%% visualize single trials

figure;
for it = 1:size(SgmLFP{1},3);
    
    subplot(7,1,1:3);
    imagesc(tx-1,fx2,squeeze(10*log10(SgmLFP2{1}(:,:,it)))');
    axis xy;ylim([30 200]);
    xlim([min(tx-1) max(tx)-1]);
    
    subplot(7,1,4:5);
    imagesc(tx-1,fx,squeeze(10*log10(SgmLFP{1}(:,:,it)))');
    axis xy;ylim([0 30]);
    xlim([min(tx-1) max(tx)-1]);
    
    subplot(7,1,6);
    selIx = find(spkDat{1}.trial{1} == trlENC(it));
    plot(spkDat{1}.time{1}(selIx),ones(1,length(spkDat{1}.timestamp{1}(selIx))),'k.');
    xlim([min(tx-1) max(tx)-1]);
    ylim([-.1 1.5]);
    axis off;
    
    subplot(7,1,7);
    hold on;
    x = dum.trial{it};
    z1 = (x-mean(x))./std(x);
    plot(bbfilt.time{1},z1,'k'); 
    x = bbfilt.trial{it};
    plot(bbfilt.time{1},(x-mean(x))./std(x),'r');     
%     for zt = 1:length(mix1{it})
%         plot([tx2(mix1{it}(zt)) tx2(mix1{it}(zt))],[-5 5],'b--');
%     end;
%     for zt = 1:length(mix2{it})
%         plot([tx2(mix2{it}(zt)) tx2(mix2{it}(zt))],[-5 5],'g--');
%     end;
%     %axis tight;
    xlim([min(tx-1) max(tx-1)]);
    ylim([-5 5]);
    axis off;
    pause;clf;
end;


%% visualize spectrogram and coherencegram

figure;
subplot(3,1,1);
a =gca;
hold on;
imagesc(tx-1,fx,10*log10(SgmLFP)');
subplot(3,1,2);
a = [a gca];
hold on;
imagesc(tx-1,fx,10*log10(SgmSPK)');
subplot(3,1,3);
a = [a gca];
hold on;
imagesc(tx-1,fx,CgmLFPSPK');

for it = 1:length(a)
    axis(a(it),'xy');
    plot(a(it),[0 0],[min(fx) max(fx)],'w');
    plot(a(it),[2 2],[min(fx) max(fx)],'w');
    axis(a(it),'tight');
    xlabel(a(it),'Time [s]');
    ylabel(a(it),'Frequency [Hz]');
    set(a(it),'ylim',[0 20]);
end;

%%
selIx = 1:length(trlENC);%find(y >= quantile(y,.6));
trlENC2 = trlENC(selIx);

for jt = 2%1:length(spkDat)
    
    cfg                     = [];
    cfg.keeptrials          = 'yes';
    cfg.channel             = spkDat{jt}.hdr.label;
    cfg.preproc.demean       = 'yes';
    cfg.preproc.detrend      = 'yes';
    cfg.trials               = selIx;
    
    [tlck] = ft_timelockanalysis( cfg ,MWdat );
    
    cfg                     =[];
    cfg.latency             = [-1 5];
    
    [tlck] = ft_selectdata( cfg, tlck );
    
    data_lfp = squeeze(tlck.trial)';
    
    %dt = [min(tlck.time) max(tlck.time)].*1e3;
    %dt = dt(1):dt(2);
    dt = [min(tlck.time) max(tlck.time)];
    dt = dt.*1e3;
    dt = dt(1):dt(2);
    
    trl = spkDat{jt}.trial{1};
    ts = spkDat{jt}.time{1}.*1e3;
    
    data_spk  = zeros(length(dt),length(trlENC2));
    data_spk2 = struct;
    data_lfp2 = zeros(length(dt),length(trlENC2));
    
    for it = 1:length(trlENC2)
        x = ts(trl == trlENC2(it));
        data_spk(:,it) = hist(x,dt);
        data_spk2(it).times = x./1e3;
        [data_lfp2(:,it)] = interpLFP(data_lfp(:,it),data_spk(:,it),[2 3]);
    end;
    
    Fs = 1/1e-3;
    TW = (size(data_lfp,1)/Fs)*1/(size(data_lfp,1)/Fs);
    params1                  = [];
    params1.tapers           = [TW , 2*TW-1];
    params1.pad              = 0;
    params1.Fs               = Fs;
    params1.err              = 0;
    params1.trialave         = 1;
    params1.fpass            = [0 30];
    
    [C,~,~,S1,S2,fx]=coherencycpb(data_lfp,data_spk,params1);
    [C2,~,~,S3,S4,fx]=coherencycpb(data_lfp2,data_spk,params1);
    
    [pval,Cemp,Cshuf] = shuffledSFC(data_lfp2,data_spk,params1,Fs);
    
    figure;
    hold on;
    plot(fx,mean(C,2));
    
    figure;
    plot(fx,mean(C2,2));
    
    figure;
    hold on;
    plot(fx,Cemp);
    pIx = find(pval <1e-4);
    plot(fx(pIx),Cemp(pIx),'b*');
    
    smp = dt./1e3;    
    D = [-2 2];
    plt = 'k';
    err = 1;    
    T = [smp(1) smp(end)];
    
    figure;
    [s,t,E] = sta(data_spk2,data_lfp2,smp,plt,[],T,D,err);
    
    delIx1 = find(sign(s(:,t==0))==1 );
    delIx2 = find(sign(s(:,t==0))==-1 );
    
    tmp1 = mean(s(delIx1,:),1);
    tmp2 = mean(s(delIx2,:),1);
    
    for zt = 1:length(delIx1)        
        s(delIx1(zt),:) = s(delIx1(zt),:)-tmp1;        
    end;           
    for zt = 1:length(delIx2)        
        s(delIx2(zt),:) = s(delIx2(zt),:)-tmp2;        
    end;
    
    
    figure;
    plot(t,mean(s,1));
    axis tight;
    
      
    figure;
    subplot(5,1,1:3);
    N = size(s,1);
    phi = hilbert(s);
    y = cos(angle(phi));
    
    ix1 = find(sign(y(:,find(t==0)))==1);
    ix2 = find(sign(y(:,find(t==0)))==-1);
    
    [~,sIdx] = sort(y(:,find(t==0)),'descend');
    imagesc(t,1:size(s,1),y(sIdx,:));
    axis xy;
    subplot(5,1,4);
    itc = phi./abs(phi);
    itc = sum(itc,1);
    itc = abs(itc)./N;
    plot(t,itc);axis tight;
    subplot(5,1,5);
    hold on;
    plot(t,mean(cos(angle(phi)),1));
    caxis([-1 1]);
    axis tight;
    
end;

%%
addpath(genpath('/home/rouxf/tbx/chronux_2_11/'));

ntrl = length(trlENC);

%toi = [min(MWdat.time{1}) max(MWdat.time{1})];
toi = [min(MWdat.time{1}) max(MWdat.time{1})];

for lt = 1:length(MWdat.label)
    
    cfg                      = [];
    cfg.channel              = MWdat.label(lt);
    cfg.keeptrials           = 'yes';
    cfg.preproc.demean       = 'yes';
    cfg.preproc.detrend      = 'yes';
    
    [dum] = ft_timelockanalysis( cfg , MWdat );
    
    cfg                     = [];
    cfg.latency             = [toi(1) toi(2)];
    
    [dum] = ft_selectdata( cfg , dum );
    
    nsamp = length(dum.time);
    tt = dum.time.*1e3;
    dt = dum.time(2)-dum.time(1);
    
    [dum] = squeeze(dum.trial)';
    
    data1 = cell(1,length(unique(uIx(:,1))));
    data2 = cell(1,length(unique(uIx(:,1))));
    cnt = 0;
    for it = 1:length(spkDat)
        for jt = 1:length(spkDat{it}.unit)
            cnt = cnt+1;
            data2{cnt} = zeros(nsamp,ntrl);
            data1{cnt} = zeros(nsamp,ntrl);
            for kt = 1:ntrl
                ts = spkDat{it}.time{jt}(spkDat{it}.trial{jt} == trlENC(kt));
                ts = ts(ts>toi(1) & ts <toi(2));
                ts = ts.*1e3;
                [n,~] = hist( ts , tt );
                data2{cnt}(:,kt) = n;
                
                [data1{cnt}(:,kt)] = interpLFP(dum(:,kt),data2{cnt}(:,kt),[3 5]);
                %data1{cnt}(:,kt) = dum(:,kt);
            end;
        end;
    end;
    
    Fs = 1/dt;
    TW = 5;%(size(data1,1)/Fs)*1/(size(data1,1)/Fs);
    params1                  = [];
    params1.tapers           = [TW , 2*TW-1];
    params1.pad              = -1;
    params1.Fs               = 1/dt;
    params1.err              = 0;
    params1.trialave         = 1;
    params1.fpass            = [0 200];
    
    movingwin = [1 .067];
    TW = movingwin(1)*1/movingwin(1);
    params2                  = [];
    params2.tapers           = [TW , 2*TW-1];
    params2.pad              = -1;
    params2.Fs               = 1/dt;
    params2.err              = 0;
    params2.trialave         = 1;
    params2.fpass            = [0 200];
    
    n = length( data2 );
    %k = 0;figure;ax1 = zeros(n*3,1);for it = 1:length( data2 );for jt=1:3;k = k+1;subplot(n,3,k);ax1(k) =gca;end;end;
    
    k = 0;
    
    for it = 1:length( data2 )
        
        [C,~,~,~,~,f3]=coherencycpb(data1{it},data2{it},params1);
        [Cgrm,~,~,~,~,t4,f4,]=cohgramcpb(data1{it},data2{it},movingwin,params2);
        
        [Sp1,f1]=mtspectrumc(data1{it},params1);
        
%         selIx = find(f1 >20 & f1 <100);
%         f1 = f1(selIx);
%         Sp1 = Sp1(selIx);
%         
%         b = regress(log10(Sp1),[log10(f1') ones(size(f1'))]);
%         yp = b(2)+b(1)*log10(f1');
%         
%         yc = log10(Sp1)-yp;
        
        [Sp2,f2,R]=mtspectrumpb(data2{it},params1);
        
        figure;
        k=1;
        subplot(131);
        %k = k+1;
        ax1(1) = gca;
        plot(ax1(k),f1,log10(Sp1));
        %plot(ax1(k),f1,yc);
        axis tight;
        %if it ==1
        title(MWdat.label(lt));
        %end;
        subplot(132);
        ax1(1) = gca;
        %k = k+1;
        plot(ax1(k),f2,Sp2./R);
        %if it ==1
        title(['Snn, Rate:',num2str(R)]);
        %end;
        axis tight;
        subplot(133);
        ax1(1) = gca;
        %k = k+1;
        plot(ax1(k),f3,C);
        axis tight;
        
        %     figure;
        %     hold on;
        %     imagesc(t-2,f2,Cgrm');axis xy;
        %     plot([0 0],[min(f2) max(f2)],'w');
        %     plot([2 2],[min(f2) max(f2)],'w');
        %     axis tight;
        %     xlabel('Time [s]');
        %     ylabel('Frequency [s]');
        %     title([chLab{it},':',num2str(sum(sum(data2{it})))]);
        
    end;
    %set(ax1(k),'XLim',[0 100]);
    
    cnt = 0;
    for jt = 1:length(spkDat)
        
        for kt = 1:length(spkDat{jt}.unit)
            
            cnt = cnt+1;
            
            data_spk = [];
            for it = 1:length(trlENC)
                selIx = find(spkDat{jt}.trial{kt} == trlENC(it));
                data_spk(it).times = spkDat{jt}.time{kt}(selIx);
            end;
            
            data_lfp = data1{cnt};
            smp = tt./1e3;
            
            D = [-.5 .5];
            plt = 'n';
            err = 1;
            
            T = [toi(1) toi(2)];
            [s,t,E] = sta(data_spk,data_lfp,smp,plt,[],T,D,err);
            %ylabel(chLab{cnt});
                        
            figure;
            subplot(6,1,1:3);
            N = size(s,1);
            phi = hilbert(s);
            imagesc(t,1:size(s,1),cos(angle(phi)));
            axis xy;
            subplot(6,1,4);
            itc = phi./abs(phi);
            itc = sum(itc,1);
            itc = abs(itc)./N;
            plot(t,itc);axis tight;
            subplot(6,1,5);
            hold on;
            plot(t,mean(cos(angle(phi)),1));
            caxis([-1 1]);
            axis tight;
            subplot(6,1,6);
            hold on;
            plot(t,mean(s,1),'k');
            plot(t,zeros(1,length(t)),'b');
            plot(t,mean(E)*ones(1,length(t)),'b');
            plot(t,-mean(E)*ones(1,length(t)),'b');
            axis tight;
            
        end;
    end;
end;



%%
rms = [];
for it = 1:length( MACdat.trial )
    
    rms(it,:) = sqrt(sum(MACdat.trial{it}.^2,2));
    
end;

%%
ms = mean(rms,1);
% 
% C = [0 0 0; .9 0 0; 0 0 .9; 0 .9 0];
% 
% figure;
idx = 1:8;
ixRef = [];
c = 0;
for it = 1:8:size(rms,2)
    c = c+1;
    
    Y = gradient(rms(idx));
    [~,mIx] = max(Y);
    
    ixRef(c) = idx(mIx);
    
    %hold on;
    %plot(idx,Y,'s-','Color',C(c,:));
    %plot(idx(mIx),Y(mIx),'ys','MarkerFaceColor','y');
    
    idx = idx+8;
    
end;
% x = MACdat.label;
% for it = 1:length(x);x{it}(regexp(x{it},' ')) = [];end;
% for it = 1:length(x);x{it}(regexp(x{it},'_')) = [];end;
% 
% ylabel('\Delta-RMS');
% set(gca,'XTick',2:2:length(MACdat.label));
% set(gca,'XTickLabel',MACdat.label(2:2:end));
% set(gca,'XTickLabelRotation',-45);

%%
cfg                     = [];
cfg.channel             = MACdat.label(ixRef);

[refDat] = ft_selectdata( cfg ,MACdat );

cfg                     = [];
cfg.avgoverchan         = 'yes';

[refDat] = ft_selectdata( cfg , refDat );

cfg                     = [];
cfg.channel             = {'*1'};

[chanDat] = ft_selectdata( cfg ,MACdat );

refMacDat = chanDat;
refMacDat.trial = {};

for it = 1:length(chanDat.trial)
    fprintf([num2str(it),'/',num2str(length(chanDat.trial))]);
    refMacDat.trial{it} = chanDat.trial{it} - repmat(refDat.trial{1},[4 1]);
    fprintf('\n');
end;
for it = 1:length(refMacDat.label)
    refMacDat.label(it) = {[chanDat.label{it},'-',refDat.label{1}]};
end;

% %%
% cfg                     = [];
% cfg.bpfilter            = 'yes';
% cfg.bpfreq              = [1 30];
% cfg.bpfilttype          = 'fir';
% 
% [filtDat] = ft_preprocessing( cfg ,refMacDat );
% 
% 
% 
%%
% [dum] = ft_appenddata( [], refDat, chanDat, refMacDat);
% 
% cfg                     = [];
% cfg.latency             = [-.1 4.5];
% 
% [dum] = ft_selectdata( cfg , dum );
% 
% cfg                     = [];
% cfg.viewmode            = 'vertical';
% 
% ft_databrowser( cfg , dum );

% %%
% cfg                     = [];
% cfg.latency             = [-1 max(refMacDat.time{1})];
% cfg.channel             = MWdat.label;
% 
% [MWdat] = ft_selectdata( cfg, MWdat );
% 
% cfg                     = [];
% cfg.bpfilter            ='yes';
% cfg.bpfreq              = [40 90];
% 
% [filt] = ft_preprocessing( cfg , MWdat );
% 
% figure;
% hold on;
% for it = 1:length( MWdat.trial )
%     
%     [mix1] = local_max(filt.trial{it});
%     [mix2] = local_max(MWdat.trial{it});
%     
%     mix = intersect(mix1,mix2);
%     
%     plot(MWdat.time{1},MWdat.trial{it},'k');
%     plot(filt.time{1},filt.trial{it},'r');
%     plot(filt.time{1}(mix),filt.trial{it}(mix),'b.');
%     xlim([2 4]);
%     
%     pause;
%     clf;
% end;

%%
cfg                     = [];
cfg.latency             = [-1 max(refMacDat.time{1})];
cfg.channel             = MWdat.label;

[MWdat] = ft_selectdata( cfg, MWdat );

dt = 250;

cfg                     = [];
cfg.bpfilter            ='yes';
cfg.bpfreq              = [3 6];
cfg.bpfilttype          = 'fir';

[filt] = ft_preprocessing( cfg , MWdat );

pIX1 = {};
epch1 = [];
c = 0;
for it = 1:length( filt.trial )
    fprintf([num2str(it),'/',num2str(length(filt.trial))]);
    
    [mix1] = local_max(filt.trial{it});
    [mix2] = local_max(MWdat.trial{it});
    mix1 =mix1;%intersect(mix1,mix2);
    
%     delIx = unique([find(diff(mix1).*1e-3 < 1/cfg.bpfreq(2))+1 find(diff(mix1).*1e-3 > 1/cfg.bpfreq(1))+1 find(sign(filt.trial{it}(mix1))~=1)]);
%     mix1(delIx) = [];
    
    pIX1{it} = mix1;
    
%     delIx = unique([find(diff(mix2).*1e-3 < 1/cfg.bpfreq(2))+1 find(diff(mix2).*1e-3 > 1/cfg.bpfreq(1))+1 find(sign(-filt.trial{it}(mix2))~=1)]);
%     mix2(delIx) = [];
    
    for jt = 1:length(mix1)
        selIx = mix1(jt)-dt:mix1(jt)+dt;
        
        if (min(selIx)>0) && (max(selIx) < length(MWdat.trial{it}))
            c = c + 1;
            epch1(c,:) = MWdat.trial{it}(selIx);
        end;
        
    end;
    
    fprintf('\n');
end;

cfg                     = [];
cfg.latency             = [-1 max(refMacDat.time{1})];
cfg.channel             = refMacDat.label(3)

[refMacDat] = ft_selectdata( cfg, refMacDat );

cfg                     = [];
cfg.bpfilter            ='yes';
cfg.bpfreq              = [3 6];
cfg.bpfilttype          = 'fir';

[filt] = ft_preprocessing( cfg , refMacDat );

pIX2 = {};
epch2 = [];
c = 0;
for it = 1:length( filt.trial )
    fprintf([num2str(it),'/',num2str(length(filt.trial))]);
    
    [mix1] = local_max(filt.trial{it});
    [mix2] = local_max(refMacDat.trial{it});
    mix1 =mix1;%intersect(mix1,mix2);
    
%     delIx = unique([find(diff(mix1).*1e-3 < 1/cfg.bpfreq(2))+1 find(diff(mix1).*1e-3 > 1/cfg.bpfreq(1))+1 find(sign(filt.trial{it}(mix1))~=1)]);
%     mix1(delIx) = [];
    
    pIX2{it} = mix1;
    
%     delIx = unique([find(diff(mix2).*1e-3 < 1/cfg.bpfreq(2))+1 find(diff(mix2).*1e-3 > 1/cfg.bpfreq(1))+1 find(sign(-filt.trial{it}(mix2))~=1)]);
%     mix2(delIx) = [];
    
    for jt = 1:length(mix1)
        selIx = mix1(jt)-dt:mix1(jt)+dt;
        
        if (min(selIx)>0) && (max(selIx) < length(refMacDat.trial{it}))
            c = c + 1;
            epch2(c,:) = refMacDat.trial{it}(selIx);
        end;
        
    end;
    
    fprintf('\n');
end;

%%
for it = 1:length( spkDat )
    
    if ~isempty(spkDat{it}.unit{1})
        for jt = 1:length( spkDat{it}.unit )
            
            x = {};
            c = 0;
            for kt = 1:length( trlENC )
                
                pt = MWdat.time{1}(pIX1{kt}).*1e3;
                
                trl = spkDat{it}.trial{jt};
                selIx = find(trl == trlENC(kt));
                
                ts = spkDat{it}.time{jt}(selIx).*1e3;
                
                for nt = 1:length(pt)
                    c =c+1;
                    dt = pt(nt) - 250:1:pt(nt) + 250;
                    selIx = find( ts >= dt(1) & ts <=dt(end));
                    x{c} = ts(selIx)-pt(nt);
                    if min(x{c})<-250
                        error('toto');
                    end;
                    if max(x{c})>250
                        error('toto');
                    end;
                end;
                
                dt = -250:250;
                h = zeros(length(x),length(dt));
                for nt = 1:length( x )
                    h(nt,:) = hist(x{nt},dt);
                end;
            end;
            
            figure;
            subplot(5,1,1:3);
            hold on;
            for nt = 1:length( x )
                spk = x{nt};
                trl = nt*ones(1,length(spk));
                spk = [spk;spk];
                trl = [trl-.5;trl+.5];
                
                line(spk,trl,'Color','k');
            end;
            xlim([min(dt) max(dt)]);
            ylim([0 length(x)+1]);
            title([num2str(it),'/',num2str(jt)]);
            subplot(5,1,4);
            bar(dt,sum(h,1));
            axis tight;
            subplot(5,1,5);
            plot(dt,mean(epch1,1));
            axis tight;
            
        end;
        
    end;
    
end;

%%
for it = 1:length( spkDat)
    
    if ~isempty(spkDat{it}.unit{1})
                
        for jt = 1:length( spkDat{it}.unit )
            LFP = [];
            c = 0;
            for kt  = 1:length( trlENC )
                
                selIx = find( spkDat{it}.trial{jt} == trlENC(kt));
                
                ts = spkDat{it}.time{jt}(selIx);
                ts = ts.*1e3;
                
                for nt = 1:length( ts )
                    ix = nearest(MWdat.time{1}.*1e3,ts(nt));
                    ix = ix-250:ix+250;
                    
                    if (min(ix)>=1) && (max(ix) <= length(MWdat.time{1}))
                        c = c+1;
                        LFP(c,:) = MWdat.trial{kt}(ix);
                    end;
                end;
                
            end;
            
            figure;
            subplot(5,1,1:4);
            imagesc(-250:250,1:size(LFP,1),LFP);
            title(num2str(length(spkDat{it}.time{jt})));
            axis xy;
            subplot(5,1,5);
            plot(-250:250,mean(LFP,1));
            
        end;
        
        
    end;
    
end;

%%
dt1 = [min(MWdat.time{1}):.1:max(MWdat.time{1})];

n1 = zeros(1,length(dt1));
for it = 1:length(pIX1)
    [x,~] = hist(MWdat.time{1}([pIX1{it}]),dt1);
    n1 = n1+x;
end;
n1 = n1./it;

dt2 = [min(MWdat.time{1}):.1:max(MWdat.time{1})];
n2 = zeros(1,length(dt2));
for it = 1:length(pIX2)
    [x,~] = hist(refMacDat.time{1}([pIX2{it}]),dt2);
    n2 = n2+x;
end;
n2 = n2./it;


figure;
subplot(121);
hold on;
bar(dt1,n1);
axis tight;
subplot(122);
hold on;
bar(dt2,n2);
axis tight;

%%

allE1 = [pIX1{:}];

t1 = MWdat.time{1}(allE1);

dt1 = min(MWdat.time{1}):.1:max(MWdat.time{1});

[n1,x1] = hist(t1,dt1);

figure;
hold on;
plot(dt1,cumsum(n1./sum(n1)),'b');
plot([0 0],[0 1],'k');
plot([2 2],[0 1],'k');
plot([min(dt1) max(dt1)],[.5 .5],'r--');
axis tight;

%%
dum1 = struct;
dum1.label = {'sta_chan_post'};
dum1.fsample = 1e3;
dum1.cfg = [];
for it = 1:size( epch1 ,1)
    dum1.trial{it} = [zeros(1,length(epch1(it,:))*9) epch1(it,:) zeros(1,length(epch1(it,:))*9)];
    dum1.time{it} = 0:1/dum1.fsample:(length(dum1.trial{it})-1)/dum1.fsample;
end;

dum2 = struct;
dum2.label = {'sta_chan_post'};
dum2.fsample = 1e3;
dum2.cfg = [];
for it = 1:size( epch2 ,1)
    dum2.trial{it} = [zeros(1,length(epch2(it,:))*9) epch2(it,:) zeros(1,length(epch2(it,:))*9)];
    dum2.time{it} = 0:1/dum2.fsample:(length(dum2.trial{it})-1)/dum2.fsample;
end;


cfg                     =[];
cfg.method              = 'mtmfft';
cfg.pad                 = 'maxperlen';
cfg.taper               = 'dpss';
cfg.tapsmofrq           = 5;%1/dum1.time{1}(end);%
cfg.keeptrials          = 'yes';

[pow1] = ft_freqanalysis(cfg,dum1);
[pow2] = ft_freqanalysis(cfg,dum2);

% cfg                     = [];
% cfg.frequency           = [3 20];%
% 
% [pow1] = ft_selectdata( cfg , pow1 );
% [pow2] = ft_selectdata( cfg , pow2 );

%%
for it = 1:size(epch1,1)
    x = epch1(it,:);    
    epch1(it,:) = (x-min(x))./(max(x)-min(x));    
end;

for it = 1:size(epch2,1)
    x = epch2(it,:);    
    epch2(it,:) = (x-min(x))./(max(x)-min(x));    
end;

%%
dum1 = struct;
dum1.label = {'sta_chan_post'};
dum1.fsample = 1e3;
dum1.cfg = [];
for it = 1:size( epch1 ,1)
    dum1.trial{it} = [zeros(1,length(epch1(it,:))*9) epch1(it,:) zeros(1,length(epch1(it,:))*9)];
    dum1.time{it} = 0:1/dum1.fsample:(length(dum1.trial{it})-1)/dum1.fsample;
end;

dum2 = struct;
dum2.label = {'sta_chan_post'};
dum2.fsample = 1e3;
dum2.cfg = [];
for it = 1:size( epch2 ,1)
    dum2.trial{it} = [zeros(1,length(epch2(it,:))*9) epch2(it,:) zeros(1,length(epch2(it,:))*9)];
    dum2.time{it} = 0:1/dum2.fsample:(length(dum2.trial{it})-1)/dum2.fsample;
end;


cfg                     =[];
cfg.method              = 'mtmfft';
cfg.pad                 = 'maxperlen';
cfg.taper               = 'dpss';
cfg.tapsmofrq           = 5;%1/dum1.time{1}(end);%
cfg.keeptrials          = 'yes';

[pow3] = ft_freqanalysis(cfg,dum1);
[pow4] = ft_freqanalysis(cfg,dum2);

% cfg                     = [];
% cfg.frequency           = [3 20];%
% 
% [pow3] = ft_selectdata( cfg , pow3 );
% [pow4] = ft_selectdata( cfg , pow4 );

%%
cfg                     = [];
cfg.frequency           = [4 9];%
cfg.avgoverfreq         = 'yes';
% 
[dum1] = ft_selectdata( cfg , pow1 );
[dum2] = ft_selectdata( cfg , pow2 );

%%
selIx1 = find(dum1.powspctrm >quantile(dum1.powspctrm,.25));
[~,sIx1] = sort(dum1.powspctrm);

selIx2 = find(dum2.powspctrm >quantile(dum2.powspctrm,.25));
[~,sIx2] = sort(dum2.powspctrm);

figure;
subplot(221);
hold on;
imagesc(-dt:dt,1:length(selIx1),epch1(selIx1,:));
x = mean(epch1(selIx1,:),1);
x = x-mean(x);
x = round(length(selIx1)/2)+round(length(selIx1)/1).*x;
plot(-dt:dt,x,'r','LineWidth',3);
axis xy;axis tight;
title(MWdat.label);
plot([-200 -200],[1 length(selIx1)],'w');
plot([200 200],[1 length(selIx1)],'w');
xlabel('Time [ms]');
ylabel('Trial #');

subplot(222);
X = log10(pow1.freq)';
%X(1) = [];
Y =squeeze(mean(log10(pow1.powspctrm(selIx1,:)),1))';
%Y(1) = [];
b = regress(Y,[X ones(size(X))]);
yp = b(2)+b(1).*X;
Yc = Y-yp;
%plot(pow1.freq,Yc,'r')
%plot(X,Y);
hold on;
%plot(X,yp,'r');
x1 = pow1.freq;
y1 = squeeze(mean((pow1.powspctrm(selIx1,:)),1));
x2 = pow3.freq;
y2 = squeeze(mean((pow3.powspctrm(selIx1,:)),1));

[ax,h1,h2] = plotyy(x1,y1,x2,y2);
axis(ax,'tight');
xlabel('Frequency [Hz]');
ylabel('Power [a.u.]');
set(ax,'XLim',[20 200]);

subplot(223);
hold on;
imagesc(-dt:dt,1:length(selIx2),epch2(selIx2,:));
x = mean(epch2(selIx2,:),1);
x = x-mean(x);
x = round(length(selIx2)/2)+round(length(selIx2)/1).*x;
plot(-dt:dt,x,'r','LineWidth',3);
axis xy;axis tight;
title(refMacDat.label);
plot([-200 -200],[1 length(selIx2)],'w');
plot([200 200],[1 length(selIx2)],'w');
xlabel('Time [ms]');
ylabel('Trial #');

subplot(224);
X = log10(pow2.freq)';
%X(1) = [];
Y =squeeze(mean(log10(pow2.powspctrm(selIx2,:)),1))';
%Y(1) = [];
b = regress(Y,[X ones(size(X))]);
yp = b(2)+b(1).*X;
Yc = Y-yp;
%plot(pow2.freq,Yc,'r');
%plot(X,Y);
%hold on;
%plot(X,yp,'r');
%plot(pow2.freq,squeeze(mean((pow2.powspctrm(selIx2,:)),1)));
x1 = pow2.freq;
y1 = squeeze(mean((pow2.powspctrm(selIx2,:)),1));
x2 = pow4.freq;
y2 = squeeze(mean((pow4.powspctrm(selIx2,:)),1));
[ax,h1,h2] = plotyy(x1,y1,x2,y2);
axis(ax,'tight');
xlabel('Frequency [Hz]');
ylabel('Power [a.u.]');
set(ax,'XLim',[20 200]);

%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
dt = 250;

cfg                     = [];
cfg.bpfilter            ='yes';
cfg.bpfreq              = [40 90];

[filt] = ft_preprocessing( cfg , MWdat );

pIX1 = {};
epch1 = [];
c = 0;
for it = 1:length( filt.trial )
    fprintf([num2str(it),'/',num2str(length(filt.trial))]);
    
    [mix1] = local_max(filt.trial{it});
    [mix2] = local_max(MWdat.trial{it});
    mix1 =intersect(mix1,mix2);
    
%     delIx = unique([find(diff(mix1).*1e-3 < 1/cfg.bpfreq(2))+1 find(diff(mix1).*1e-3 > 1/cfg.bpfreq(1))+1 find(sign(filt.trial{it}(mix1))~=1)]);
%     mix1(delIx) = [];
    
    pIX1{it} = mix1;
    
%     delIx = unique([find(diff(mix2).*1e-3 < 1/cfg.bpfreq(2))+1 find(diff(mix2).*1e-3 > 1/cfg.bpfreq(1))+1 find(sign(-filt.trial{it}(mix2))~=1)]);
%     mix2(delIx) = [];
    
    for jt = 1:length(mix1)
        selIx = mix1(jt)-dt:mix1(jt)+dt;
        
        if (min(selIx)>0) && (max(selIx) < length(MWdat.trial{it}))
            c = c + 1;
            epch1(c,:) = MWdat.trial{it}(selIx);
        end;
        
    end;
    
    fprintf('\n');
end;


cfg                     = [];
cfg.bpfilter            ='yes';
cfg.bpfreq              = [40 90];

[filt] = ft_preprocessing( cfg , refMacDat );

pIX2 = {};
epch2 = [];
c = 0;
for it = 1:length( filt.trial )
    fprintf([num2str(it),'/',num2str(length(filt.trial))]);
    
    [mix1] = local_max(filt.trial{it});
    [mix2] = local_max(refMacDat.trial{it});
    mix1 =intersect(mix1,mix2);
    
%     delIx = unique([find(diff(mix1).*1e-3 < 1/cfg.bpfreq(2))+1 find(diff(mix1).*1e-3 > 1/cfg.bpfreq(1))+1 find(sign(filt.trial{it}(mix1))~=1)]);
%     mix1(delIx) = [];
    
    pIX2{it} = mix1;
    
%     delIx = unique([find(diff(mix2).*1e-3 < 1/cfg.bpfreq(2))+1 find(diff(mix2).*1e-3 > 1/cfg.bpfreq(1))+1 find(sign(-filt.trial{it}(mix2))~=1)]);
%     mix2(delIx) = [];
    
    for jt = 1:length(mix1)
        selIx = mix1(jt)-dt:mix1(jt)+dt;
        
        if (min(selIx)>0) && (max(selIx) < length(refMacDat.trial{it}))
            c = c + 1;
            epch2(c,:) = refMacDat.trial{it}(selIx);
        end;
        
    end;
    
    fprintf('\n');
end;

%%
dum1 = struct;
dum1.label = {'sta_chan_post'};
dum1.fsample = 1e3;
dum1.cfg = [];
for it = 1:size( epch1 ,1)
    dum1.trial{it} = [zeros(1,length(epch1(it,:))*1) epch1(it,:) zeros(1,length(epch1(it,:))*1)];
    dum1.time{it} = 0:1/dum1.fsample:(length(dum1.trial{it})-1)/dum1.fsample;
end;

dum2 = struct;
dum2.label = {'sta_chan_post'};
dum2.fsample = 1e3;
dum2.cfg = [];
for it = 1:size( epch2 ,1)
    dum2.trial{it} = [zeros(1,length(epch2(it,:))*1) epch2(it,:) zeros(1,length(epch2(it,:))*1)];
    dum2.time{it} = 0:1/dum2.fsample:(length(dum2.trial{it})-1)/dum2.fsample;
end;


cfg                     =[];
cfg.method              = 'mtmfft';
cfg.pad                 = 'maxperlen';
cfg.taper               = 'dpss';
cfg.tapsmofrq           = 5;%1/dum1.time{1}(end);%
cfg.keeptrials          = 'yes';

[pow1] = ft_freqanalysis(cfg,dum1);
[pow2] = ft_freqanalysis(cfg,dum2);

% cfg                     = [];
% cfg.frequency           = [20 200];%[3 20];%
% 
% [pow1] = ft_selectdata( cfg , pow1 );
% [pow2] = ft_selectdata( cfg , pow2 );

%%
for it = 1:size(epch1,1)
    x = epch1(it,:);    
    epch1(it,:) = (x-min(x))./(max(x)-min(x));    
end;

for it = 1:size(epch2,1)
    x = epch2(it,:);    
    epch2(it,:) = (x-min(x))./(max(x)-min(x));    
end;

%%
dum1 = struct;
dum1.label = {'sta_chan_post'};
dum1.fsample = 1e3;
dum1.cfg = [];
for it = 1:size( epch1 ,1)
    dum1.trial{it} = [zeros(1,length(epch1(it,:))*1) epch1(it,:) zeros(1,length(epch1(it,:))*1)];
    dum1.time{it} = 0:1/dum1.fsample:(length(dum1.trial{it})-1)/dum1.fsample;
end;

dum2 = struct;
dum2.label = {'sta_chan_post'};
dum2.fsample = 1e3;
dum2.cfg = [];
for it = 1:size( epch2 ,1)
    dum2.trial{it} = [zeros(1,length(epch2(it,:))*1) epch2(it,:) zeros(1,length(epch2(it,:))*1)];
    dum2.time{it} = 0:1/dum2.fsample:(length(dum2.trial{it})-1)/dum2.fsample;
end;


cfg                     =[];
cfg.method              = 'mtmfft';
cfg.pad                 = 'maxperlen';
cfg.taper               = 'dpss';
cfg.tapsmofrq           = 5;%1/dum1.time{1}(end);%
cfg.keeptrials          = 'yes';

[pow3] = ft_freqanalysis(cfg,dum1);
[pow4] = ft_freqanalysis(cfg,dum2);
% 
% cfg                     = [];
% cfg.frequency           = [20 200];%[3 20];%
% 
% [pow3] = ft_selectdata( cfg , pow3 );
% [pow4] = ft_selectdata( cfg , pow4 );

%%
cfg                     = [];
cfg.frequency           = [40 90];%[3 20];%
cfg.avgoverfreq         = 'yes';

[dum1] = ft_selectdata( cfg , pow1 );
[dum2] = ft_selectdata( cfg , pow2 );

%%
selIx1 = find(dum1.powspctrm >quantile(dum1.powspctrm,.995));
[~,sIx1] = sort(dum1.powspctrm);

selIx2 = find(dum2.powspctrm >quantile(dum2.powspctrm,.995));
[~,sIx2] = sort(dum2.powspctrm);

figure;
subplot(221);
hold on;
imagesc(-dt:dt,1:length(selIx1),epch1(selIx1,:));
x = mean(epch1(selIx1,:),1);
x = x-mean(x);
x = round(length(selIx1)/2)+round(length(selIx1)/1).*x;
plot(-dt:dt,x,'r','LineWidth',3);
axis xy;axis tight;
title(MWdat.label);
plot([-20 -20],[1 length(selIx1)],'w');
plot([20 20],[1 length(selIx1)],'w');
xlabel('Time [ms]');
ylabel('Trial #');
xlim([-250 250]);

subplot(222);
X = log10(pow1.freq)';
%X(1) = [];
Y =squeeze(mean(log10(pow1.powspctrm(selIx1,:)),1))';
%Y(1) = [];
b = regress(Y,[X ones(size(X))]);
yp = b(2)+b(1).*X;
Yc = Y-yp;
%plot(pow1.freq,Yc,'r')
%plot(X,Y);
hold on;
%plot(X,yp,'r');
x1 = pow1.freq;
y1 = squeeze(mean((pow1.powspctrm(selIx1,:)),1));
x2 = pow3.freq;
y2 = squeeze(mean((pow3.powspctrm(selIx1,:)),1));

[ax,h1,h2] = plotyy(x1,y1,x2,y2);
axis(ax,'tight');
xlabel('Frequency [Hz]');
ylabel('Power [a.u.]');
set(ax,'XLim',[0 20]);

subplot(223);
hold on;
imagesc(-dt:dt,1:length(selIx2),epch2(selIx2,:));
x = mean(epch2(selIx2,:),1);
x = x-mean(x);
x = round(length(selIx2)/2)+round(length(selIx2)/1).*x;
plot(-dt:dt,x,'r','LineWidth',3);
axis xy;axis tight;
title(refMacDat.label);
plot([-20 -20],[1 length(selIx2)],'w');
plot([20 20],[1 length(selIx2)],'w');
xlabel('Time [ms]');
ylabel('Trial #');
xlim([-250 250]);

subplot(224);
X = log10(pow2.freq)';
%X(1) = [];
Y =squeeze(mean(log10(pow2.powspctrm(selIx2,:)),1))';
%Y(1) = [];
b = regress(Y,[X ones(size(X))]);
yp = b(2)+b(1).*X;
Yc = Y-yp;
%plot(pow2.freq,Yc,'r');
%plot(X,Y);
%hold on;
%plot(X,yp,'r');
%plot(pow2.freq,squeeze(mean((pow2.powspctrm(selIx2,:)),1)));
x1 = pow2.freq;
y1 = squeeze(mean((pow2.powspctrm(selIx2,:)),1));
x2 = pow4.freq;
y2 = squeeze(mean((pow4.powspctrm(selIx2,:)),1));
[ax,h1,h2] = plotyy(x1,y1,x2,y2);
axis(ax,'tight');
xlabel('Frequency [Hz]');
ylabel('Power [a.u.]');
set(ax,'XLim',[0 20]);

%%

[erpDat,TFR] = compute_ERP_and_TFR_macroData_EM(refMacDat);

%%
cfg                     = [];
cfg.channel             = refMacDat.label(3);

[dum] = ft_selectdata( cfg , TFR );

dum.label = {'R PH_1'};

cfg                     = [];
cfg.frequency           = [60 70];
cfg.avgoverfreq         = 'yes';

[dum] = ft_selectdata( cfg , dum );

cfg                     = [];
cfg.demean              = 'yes';
cfg.detrend             = 'yes';
cfg.derivative          = 'yes';

[dum] = ft_preprocessing( cfg,dum );

cfg                     =[];
cfg.method              = 'mtmfft';
cfg.pad                 = 'maxperlen';
cfg.taper               = 'dpss';
cfg.tapsmofrq           = 1/dum.time{1}(end);

[pow] = ft_freqanalysis(cfg,dum);

cfg                     = [];
cfg.frequency           = [0 20];

[pow] = ft_selectdata( cfg , pow );

figure;
plot(pow.freq,pow.powspctrm);

%%
cfg                     = [];
cfg.output              = 'pow';
cfg.keeptrials          = 'yes';
cfg.toi                 = [min(TFR.time) max(TFR.time)];

[TFR,slp]=sh_subtr1of(cfg,TFR);

%%
cfg                     = [];
cfg.frequency           = [0 19];

[TFRlow] = ft_selectdata( cfg , TFR );

cfg                     = [];
cfg.frequency           = [20 max(TFR.freq)];

[TFRhigh] = ft_selectdata( cfg , TFR );

%%
ntrl  = size(TFRhigh.powspctrm,1);
nchan = size(TFRhigh.powspctrm,2);
nfreq = size(TFRhigh.powspctrm,3);
ntpt  = size(TFRhigh.powspctrm,4);

selIx = find(TFRhigh.time >=-1 & TFRhigh.time <=-.05);
M = squeeze(mean(mean(TFRhigh.powspctrm(:,:,:,selIx),4),1));
SD = squeeze(std(mean(TFRhigh.powspctrm(:,:,:,selIx),4),0,1));

Z = zeros(size(TFRhigh.powspctrm));
for it = 1:size(M,1)
    m = repmat(M(it,:),[ntpt 1])';
    sd = repmat(SD(it,:),[ntpt 1])';
    for jt = 1:size(TFRhigh.powspctrm,1)
        y = squeeze(TFRhigh.powspctrm(jt,it,:,:));
        Z(jt,it,:,:) = (y - m)./sd;
    end;
end;

%%
ntrl  = size(TFRlow.powspctrm,1);
nchan = size(TFRlow.powspctrm,2);
nfreq = size(TFRlow.powspctrm,3);
ntpt  = size(TFRlow.powspctrm,4);

selIx = find(TFRlow.time >=-1 & TFRlow.time <=-.05);
M = squeeze(mean(mean(TFRlow.powspctrm(:,:,:,selIx),4),1));
SD = squeeze(std(mean(TFRlow.powspctrm(:,:,:,selIx),4),0,1));

Z2 = zeros(size(TFRlow.powspctrm));
for it = 1:size(M,1)
    m = repmat(M(it,:),[ntpt 1])';
    sd = repmat(SD(it,:),[ntpt 1])';
    for jt = 1:size(TFRlow.powspctrm,1)
        y = squeeze(TFRlow.powspctrm(jt,it,:,:));
        Z2(jt,it,:,:) = (y - m)./sd;
    end;
end;

%%
a = [];
b = [];
ca = [];
cb = [];
for it = 1:length(TFR.label)
    figure;
    subplot(5,1,1:3);
    a(it) = gca;
    hold on;   
    
    imagesc(TFRlow.time,TFRhigh.freq,squeeze(mean((TFRhigh.powspctrm(:,it,:,:)),1)));
    axis xy;
    plot([0 0],[min(TFRhigh.freq) max(TFRhigh.freq)],'w');
    plot([2 2],[min(TFRhigh.freq) max(TFRhigh.freq)],'w');
    axis tight;
    chName = TFRhigh.label(it);
    chName{:}(regexp(chName{:},'_')) = [];
    ca(it,:) = caxis;
    title([chName{:},' (',num2str(size(TFRhigh.powspctrm,1)),' trials)']);
    
    subplot(5,1,4:5);
    b(it) = gca;
    hold on;
    imagesc(TFRlow.time,TFRlow.freq,squeeze(mean((TFRlow.powspctrm(:,it,:,:)),1)));
    axis xy;
    plot([0 0],[min(TFRlow.freq) max(TFRlow.freq)],'w');
    plot([2 2],[min(TFRlow.freq) max(TFRlow.freq)],'w');
    axis tight;
    cb(it,:) = caxis;
    
    %%
    cfg                     = [];
    cfg.avgovertime         = 'yes';
    cfg.nanmean             = 'yes';
    cfg.channel             = TFRlow.label(it);
    
    [dum1]  = ft_selectdata( cfg , TFRlow );
    
    cfg                     = [];
    cfg.avgovertime         = 'yes';
    cfg.nanmean             = 'yes';
    cfg.channel             = TFRhigh.label(it);
    
    [dum2]  = ft_selectdata( cfg , TFRhigh );
    
    cfg                     = [];
    cfg.variance            = 'yes';
    cfg.jackknife           = 'yes';
    
    [dum1]  = ft_freqdescriptives( cfg , dum1 );
    [dum2]  = ft_freqdescriptives( cfg , dum2 );
    
    %%
    figure;
    subplot(121);
    hold on;
    plot(dum1.freq,dum1.powspctrm);
    plot(dum1.freq,dum1.powspctrm-dum1.powspctrmsem);
    plot(dum1.freq,dum1.powspctrm+dum1.powspctrmsem);
    axis tight;
    title('Low frequencies');
    
    subplot(122);
    hold on;
    plot(dum2.freq,dum2.powspctrm);
    plot(dum2.freq,dum2.powspctrm-dum2.powspctrmsem);
    plot(dum2.freq,dum2.powspctrm+dum2.powspctrmsem);
    axis tight;
    title('High frequencies');
        
    
end;
for it = 1:length(a);caxis(a(it),ca(end,:));end;
for it = 1:length(b);caxis(b(it),cb(end,:));end;

%%
a = [];
b = [];
ca = [];
cb = [];
for it = 1:length(TFR.label)
    figure;
    subplot(5,1,1:3);
    a(it) = gca;
    hold on;
    imagesc(TFRlow.time,TFRhigh.freq,squeeze(mean((Z(:,it,:,:)),1)));
    axis xy;
    plot([0 0],[min(TFRhigh.freq) max(TFRhigh.freq)],'w');
    plot([2 2],[min(TFRhigh.freq) max(TFRhigh.freq)],'w');
    axis tight;
    chName = TFRhigh.label(it);
    chName{:}(regexp(chName{:},'_')) = [];
    ca(it,:) = caxis;
    title([chName{:},' (',num2str(size(TFRhigh.powspctrm,1)),' trials)']);
    
    subplot(5,1,4:5);
    b(it) = gca;
    hold on;
    imagesc(TFRlow.time,TFRlow.freq,squeeze(mean((Z2(:,it,:,:)),1)));
    axis xy;
    plot([0 0],[min(TFRlow.freq) max(TFRlow.freq)],'w');
    plot([2 2],[min(TFRlow.freq) max(TFRlow.freq)],'w');
    axis tight;
    cb(it,:) = caxis;
end;
for it = 1:length(a);caxis(a(it),ca(end,:));end;
for it = 1:length(b);caxis(b(it),cb(end,:));end;