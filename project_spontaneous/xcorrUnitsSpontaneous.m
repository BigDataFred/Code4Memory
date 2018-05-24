%%
addpath('/home/rouxf/tbx/releaseDec2015/');
addpath(genpath('/home/rouxf/tbx/wave_clus/'));
addpath('/home/rouxf/prj/Bham/code/mcode/utils/spikes/');
addpath(genpath('/home/rouxf/tbx/osort-v3-rel/'));
addpath('/home/rouxf/prj/Bham/code/mcode/params/');

%% parameters to read the CSC data
FieldSelection2 = [];
FieldSelection2(1) = 1;%timestamps
FieldSelection2(2) = 0;
FieldSelection2(3) = 0;%sample freq
FieldSelection2(4) = 0;
FieldSelection2(5) = 1;%samples

%% set data parameters
Fs =32e3;
timeStampsPerSec = (1/Fs)*1e6;

[par] = set_parameters_Bham(Fs);

[b,a] = butter(4,[15]./(1e3/2),'low');% apply low-pass for LFP

%%
pID = {'P07','P09'};
expMode = {'fVSpEM','Spontaneous'};

curPat = 2;
curExp = 2;

p2d = ['/media/rouxf/rds-share/Archive/MICRO/',pID{curPat},'/',expMode{curExp},'/'];
sesh = dir(p2d);
sesh(1:2) = [];
sesh = sesh([sesh(:).isdir]);

sesh = {sesh(:).name}

%%
%fn = {'CSC_midHippL5.ncs','CSC_midHippL6.ncs'};
fn = {'CSC_paraHippL5.ncs','CSC_paraHippL6.ncs'};

[ ~,chanN1, ~] = fileparts(fn{1});
[ ~,chanN2, ~] = fileparts(fn{2});
chanN1(regexp(chanN1,'CSC_'):regexp(chanN1,'CSC_')+3) = [];
chanN2(regexp(chanN2,'CSC_'):regexp(chanN2,'CSC_')+3) = [];

%%
for curSesh = 2%1:length(sesh)%1%
    
sortedSpikes = {};
wltCoeffs = {};
for curMW = 1:length(fn)
    
    [~, dataSamples,hdr] = Nlx2MatCSC_v3([p2d,sesh{curSesh},'/',fn{curMW}], FieldSelection2, 1, 1, []);
    
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
    [dataSamples] = double(dataSamples(:))';
    [dataSamples] = dataSamples.*scalef.*1e6;
    
    %% spike detection (2nd pass)
    par2 = par;
    par2.fnamespc = [par2.fnamespc,num2str(curMW)];
    par2.fname_in = [par2.fname_in,num2str(curMW)];
    par2.segments = 6;%floor(length(dataSamples)/(5*Fs));
    par2.interpolation = 'n';
    
    nsmp = floor(length(dataSamples)/(par2.segments));%Fs;
    ix = 1:nsmp;
    noiseTraces = [];
    spikeWaveforms = [];
    spikeTimestamps = [];
    noiseSDTmp = [];
    
    for kt = 1:par2.segments
        
        [~,spikeWaveformsTmp, ~, spikeTimestampsTmp,~,~,noiseTracesTmp] = amp_detect(dataSamples(ix),par2);
        
        spikeWaveforms  = [spikeWaveforms;spikeWaveformsTmp];
        spikeTimestamps = [spikeTimestamps ix(spikeTimestampsTmp)];
        noiseSDTmp      = [noiseSDTmp mean(std(noiseTracesTmp,0,2),1)];
        
        if size(noiseTracesTmp,1)>1
            noiseTracesTmp2=[];
            noiseTracesTmp2(1:size(noiseTracesTmp,1) ,1) = ones( size(noiseTracesTmp,1), 1 )*1i;
            noiseTracesTmp2(1:size(noiseTracesTmp,1),2:size(noiseTracesTmp,2)+1)=noiseTracesTmp;
            noiseTraces = [noiseTraces; noiseTracesTmp2];
            noiseTracesTmp2 = [];
            noiseTracesTmp = [];
        end;
        if kt < par2.segments-1
            ix = ix+nsmp;
        else
            ix = ix(end)+1:length(dataSamples);
        end;
    end;    
    
    %% decorrelate and upsample spike-waveforms
    if ( ~isempty(spikeWaveforms) ) && ( size(spikeWaveforms,1)>1 ) && ( size(spikeWaveforms,2) >1 )
        
        [trans] = posthocWhiten(noiseTraces, spikeWaveforms,[]);
        [trans] = [trans(:,1)*ones(1,2) trans trans(:,end)*ones(1,2)];
        
        par2.interpolation = 'y';
        [trans] = int_spikes(trans,par);
        
        [spikeWaveforms] = [spikeWaveforms(:,1)*ones(1,2) spikeWaveforms spikeWaveforms(:,end)*ones(1,2)];
        [spikeWaveforms] = int_spikes(spikeWaveforms,par);
        
    end;
    
    %%
    waveclus                            = [];
    waveclus.spikes                     = spikeWaveforms;
    waveclus.index                      = spikeTimestamps;
    
    [dim] = size(spikeWaveforms);
    
    %% do spike sorting
    par2.filename = [fn{curMW}];
    
    dum                     = [];
    dum.newSpikeTimes       = [];
    dum.assignedCluster     = [];
    dum.wavf                = [];
    dum.num_clus            = [];
    
    dum2                    = [];
    
    %%
    if dim(1) > dim(2)% the number of spike events must larger than the number of samples in the waveform
        [dum,dum2] = doSpikeSorting_waveclus( waveclus , par2 );
        
        dum.wavf_decorr         = trans;%spikeWaveforms2;%trans;%
        dum.SD                  = mean(noiseSDTmp(noiseSDTmp>0));
        trans = [];        
        
        %% delete noise clusters
        [delIx] = find(dum.assignedCluster==0);
        
        dum.assignedCluster( delIx ) = [];
        dum.newSpikeTimes( delIx ) = [];
        dum.wavf(delIx,:) = [];
        dum.wavf_decorr(delIx,:) = [];
        dum2(delIx,:) = [];
    end;
    
    %% save the clustered data
    sortedSpikes{curMW} = dum;
    wltCoeffs{curMW} = dum2;
end;

%%

[spkDat1] = sortedSpikes{(1)}

[spkDat2] = sortedSpikes{(2)}

cID1 = unique(spkDat1.assignedCluster)
cID2 = unique(spkDat2.assignedCluster)

p = [];
cnt = 0;
for it = 1:length(cID1)
    for jt = 1:length(cID2)
        cnt = cnt + 1;
        p(cnt,:) = [cID1(it) cID2(jt)];
    end
end;

p = [1 1];

for it = 1:size(p,1)
    
    ts1 = spkDat1.newSpikeTimes(spkDat1.assignedCluster ==p(it,1));
    ts2 = spkDat2.newSpikeTimes(spkDat2.assignedCluster ==p(it,2));
    
    wvf1 = spkDat1.wavf(spkDat1.assignedCluster ==p(it,1),:);
    wvf2 = spkDat2.wavf(spkDat2.assignedCluster ==p(it,2),:);
    
    %%
    [lfpTime] = [0:length(dataSamples)-1]./Fs;
    [lfpTime] = lfpTime.*1e3;
    ts1 = lfpTime(ts1);
    ts2 = lfpTime(ts2);
    
    %%
    dt = lfpTime(1):1:lfpTime(end);    
    [n1,~] = hist(ts1,dt);
    [n2,~] = hist(ts2,dt);
    [n3,~] = hist(sort([ts1 ts2]),dt);
    
    [xc,lag] = xcorr(n1,n1,500);
    xc(501) = 0;
    [xc2,lag] = xcorr(n3,500);
    xc2(501) = 0;
    filtXC = filtfilt(b,a,xc);
    filtXC2 = filtfilt(b,a,xc2);
    
    %%
    figure;
    subplot(221);
    hold on;
    plot(linspace(0,2,64),wvf1,'g');
    plot(linspace(0,2,64),wvf2,'Color',[.75 .75 .75]);
    %plot(linspace(0,2,64),mean(wvf1,1),'k','LineWidth',3);
    title([chanN1,'-clu',num2str(p(it,1))]);
    subplot(222);
    [dx] = diff(ts1);
    [n,~] = hist(dx,0:251);
    bar(0:250,n(1:end-1));xlim([-10 250]);
    subplot(223);
    hold on;
    plot(linspace(0,2,64),wvf2,'r');
    plot(linspace(0,2,64),mean(wvf2,1),'k','LineWidth',3);
    title([chanN2,'-clu',num2str(p(it,2))]);
    subplot(224);
    [dx] = diff(ts2);
    [n,~] = hist(dx,0:251);
    bar(0:250,n(1:end-1));
    xlim([-10 250]);
                
    figure;
    hold on;
    bar(lag,xc,'FaceColor',[.75 .75 .75],'EdgeColor',[.75 .75 .75]);  
    y = xc;
    x = lag;    
    plot(lag,filtXC,'k','LineWidth',3);   
    axis tight;
    x = [pID{curPat},'_','SE',num2str(curSesh),'-',sesh{curSesh},'_',expMode{curExp}];
    x(regexp(x,'_')) = '-';
    title(x);
    
    figure;
    hold on;
    bar(lag,xc2,'FaceColor',[.75 .75 .75],'EdgeColor',[.75 .75 .75]);        
    plot(lag,filtXC2,'k','LineWidth',3);   
    axis tight;
    x = [pID{curPat},'_','SE',num2str(curSesh),'-',sesh{curSesh},'_',expMode{curExp}];
    x(regexp(x,'_')) = '-';
    title(x);
       
    n = size(wltCoeffs{1},2);
    p = zeros(n*(n-1)/2,2);
    cnt = 0;
    for kt = 1:n
        for lt = kt+1:n
            cnt = cnt+1;
            p(cnt,:) = [kt lt];
        end;
    end;
%     
%     figure;
%     for kt = 1:size(p,1)
%         subplot(5,9,kt);
%         hold on;
%         plot(wltCoeffs{1}(:,p(kt,1)),wltCoeffs{1}(:,p(kt,2)),'bo');
%         plot(wltCoeffs{2}(:,p(kt,1)),wltCoeffs{2}(:,p(kt,2)),'rx');
%     end;
    
    dum = gcf;
    for kt = 1:dum.Number
        set(kt,'Color','w')
    end;
    
end;
end;
savepath = '/home/rouxf/tmp/figuresLNM/';
for it = 1:get(gcf,'Number')
    set(it,'Color','w');
    saveas(it,[savepath,'crossCorrelation',num2str(it),'.fig']);
end;
