%% set the path def
restoredefaultpath;
addpath(genpath('C:\toolbox\AnalysisFred\'));
addpath(genpath('C:\toolbox\Nlx\MatlabImportExport_v6.0.0\'));
addpath(genpath('C:\toolbox\osort-v3.0\osort-v3-rel\code\'));
addpath(genpath('C:\toolbox\chronux_2_12\'));
addpath('C:\toolbox\fieldtrip-20160309\');
ft_defaults;
%% open pool of workers
%parpool('SpmdEnabled',false);
clear;
clc;

%% set the data-paths directions
p2Nlxdata = 'C:\Data\Nlx\tuning\';
p2logdat = 'C:\Experiments\EMpairs_v4_2016-1007\log\Tuning\';
p2logparams = 'C:\Experiments\EMpairs_v4_2016-1007\params\';

%% select session to analyze
sesh = 3;
pre = 0.5;
post = 1.5;
chanSel = [ 3 ];

%% create the session labels of the Nlx data
Nlxdat = {};
Nlxdat{1} = 'P02_TS01_2016-Jul-10_10-50-00\';
Nlxdat{2} = 'P02_TS07_2016-Jul-18_11_56_08\';
Nlxdat{3} = 'P02_TS08_2016-Jul-19_10_54_39\';
%% create the session labels of the log data
Tlogfile = {};
Tlogfile{1} = 'P02_TS02_log_ctune_10-Jul-2016_10_46_xx.txt';
Tlogfile{2} = 'P02_TS07__log_ctune_18-Jul-2016_11_56_08.txt';
Tlogfile{3} = 'P02_TS09__log_ctune_19-Jul-2016_10_54_39.txt';

%% get the logfile data
f.p2logf = 'C:\Experiments\EMpairs_v4_2016-1007\log\Tuning\';
f.logf = Tlogfile{sesh};
p.ncols = 8;

[LogDat] = getNewLogDataTune(f,p);

%% get the event file data
ev_filename = 'Events.nev';

[event_dat.events,~] = getRawTTLs([p2Nlxdata,Nlxdat{sesh},ev_filename]);

[event_dat.ttl_idx] = find(event_dat.events(:,2)==7);

event_dat.n = length(event_dat.ttl_idx);

if LogDat.n ~= event_dat.n
    error('number of events is out of range');
end;
%% extract Fs and ADV
csc_filename = 'CSC_LA1.ncs';

[csc_dat.timestamps,~,csc_dat.hdr] = getRawCSCData( [p2Nlxdata,Nlxdat{sesh},csc_filename], [], [], 1 );
csc_dat.FirstTimeStamp = csc_dat.timestamps(1);

[idx] = getNlxHeaderParam(csc_dat.hdr,'SamplingFrequency');
csc_dat.Fs = str2double(csc_dat.hdr{idx}(20:end));

[idx] = getNlxHeaderParam(csc_dat.hdr,'ADBitVolts');
csc_dat.ADV = str2double(csc_dat.hdr{idx}(12:end))*1e6;

csc_dat.timeStampsPerSample = 1/csc_dat.Fs*1e6;

%% convert event ts to samples
event_dat.samples = zeros(size(event_dat.events,1),1);
for it = 1:size(event_dat.events,1)
    csc_dat.timestamps;
    event_dat.samples(it) = (event_dat.events(it,1)-double(csc_dat.FirstTimeStamp))./csc_dat.timeStampsPerSample + 1;
end;

if unique(event_dat.events(event_dat.ttl_idx,2)) ~=7
    error('trigger values are out of range');
end;

if length(event_dat.samples(event_dat.ttl_idx)) ~=LogDat.n
    error('number of events out of range');
end;

%%
event_dat.trl = zeros(length(event_dat.samples),3);
for it =1:length(event_dat.samples)
    event_dat.trl(it,:) = [event_dat.samples(it)-(2*pre*csc_dat.Fs) event_dat.samples(it)+(2*post*csc_dat.Fs) -(2*pre*csc_dat.Fs)];
end;
event_dat.trl = floor(event_dat.trl);
event_dat.trl(sign(event_dat.trl(:,1))==-1,:) = [];
event_dat.trl(event_dat.trl(:,2)>length(csc_dat.timestamps)*512,:) = [];

%%
files = dir([p2Nlxdata,Nlxdat{sesh},'*.ncs']);

for ct = 1%:length(chanSel)      
    
    %% st = tic;
    fprintf([num2str(ct),'/',num2str(length(chanSel))]);
    fprintf('\n');
    
    csc_filename = files(chanSel(ct)).name;
    
    [csc_dat.timestamps,csc_dat.dataSamples,~] = getRawCSCData( [p2Nlxdata,Nlxdat{sesh},csc_filename], [], [], 1 );
        
    csc_dat.dataSamples = csc_dat.dataSamples*csc_dat.ADV;
    [csc_dat.tsi] = timeStampinter(csc_dat.timestamps);
    %%
    nsamp = diff(event_dat.trl(1,1:2),[],2)+1;
    csc_dat.segmented = zeros(size(event_dat.trl,1),nsamp);
    for it = 1:size(event_dat.trl,1)        
        csc_dat.segmented(it,:) = csc_dat.dataSamples(event_dat.trl(it,1):event_dat.trl(it,2));
    end;
    %% make the LFP data
    dsf =32000/1280;
    dsFs = csc_dat.Fs/dsf;
    
    Wp = 300./csc_dat.Fs/2;
    [b,a] = butter(4,Wp);
    
    dum = csc_dat.segmented;    
    nsamp = size(dum,2);
    lfp_dat = zeros(size(dum));
    for it = 1:size(dum,1)
        
        padded = [fliplr(dum(it,:)) dum(it,:) fliplr(dum(it,:)) ];        
        padded = filtfilt(b,a,padded);
        
        lfp_dat(it,:) = padded(nsamp+1:nsamp+nsamp);
        clear padded;
    end;
    %%
    ds_lfp = zeros(size(lfp_dat,1),floor(size(lfp_dat,2)/dsf)+1);
    for it = 1:size(lfp_dat,1)
        ds_lfp(it,:) = downsample(lfp_dat(it,:),dsf);
    end;
    lfp_dat = ds_lfp';
    clear ds_lfp;
    %%
    lfp_dat = lfp_dat - repmat(mean(lfp_dat,1),[size(lfp_dat,1) 1]);
    
    data=locdetrend(lfp_dat,dsFs,[.1 .05]);
    %%
    ktapers = 8;
    NW = (ktapers+1)/2;
    
    params = [];
    params.tapers =  [NW ktapers];
    params.Fs = dsFs;
    params.fpass = [0 dsFs/2];
    params.pad = 8;
    
    for it = 1:size(lfp_dat,2)
        fprintf([num2str(it),'/',num2str(size(lfp_dat,2))]);
        [lfp_dat(:,it)]=rmlinesc(lfp_dat(:,it)',params,.05,'n');
        fprintf('\n');
    end;
    %%
    ktapers = 1;
    NW = (ktapers+1)/2;
    
    movingwin = [1 .25];
    
    params = [];
    params.Fs =  dsFs;
    params.pad = 4;
    params.tapers = [NW ktapers];
    params.fpass = [2 20];
    params.err = 0;
    params.trialave = 0;
    
    [S1,t1,f1]=mtspecgramc(lfp_dat,movingwin,params);
    %%
    ktapers = 8;
    NW = (ktapers+1)/2;
    
    movingwin = [.25 .066];
    
    params = [];
    params.Fs = dsFs;
    params.pad = 4;
    params.tapers = [NW ktapers];
    params.fpass = [20 200];
    params.err = 0;
    params.trialave = 0;
    
    [S2,t2,f2]=mtspecgramc(lfp_dat,movingwin,params);
    %%
    figure;
    subplot(7,1,1:3);
    imagesc(t2-2*pre,f2,squeeze(mean(20*log10(S2),3))');
    axis xy;
    subplot(7,1,4:6);
    imagesc(t1-2*pre,f1,squeeze(mean(20*log10(S1),3))');
    axis xy;
    subplot(717);
    plot(-2:1/dsFs:2,mean(lfp_dat,2))
    %%
    base_idx = find(t2-2*pre <=0);
    b = mean(S2(base_idx,:,:),1);
    b = median(b,3);
    b = repmat(b,[size(S2,1) 1 size(S2,3)]);
    %%
    event_dat.ttl_idx1 = find( event_dat.events(:,2)==7 );
    
    k = 0;    
    event_dat.ttl_idx2 = [];
    for it = 1:length(event_dat.ttl_idx1)
        for jt = 1:2
            
            ix = event_dat.ttl_idx1(it)+jt;
            if ( event_dat.events(ix,2)==0 )
                k = k +1;
                event_dat.ttl_idx2(k) = ix;
                break;
            end;
        end;
    end;
    
    if length(event_dat.ttl_idx1) ~= length(event_dat.ttl_idx2)
        error('number of events is out of range');
    end;
    
    params = [];
    params.Fs = csc_dat.Fs;   
    dt = 3;
    
    [event_dat.ttlavg1] = eta(dt,event_dat.samples(event_dat.ttl_idx1),csc_dat.dataSamples,params);
    [event_dat.ttlavg2] = eta(dt,event_dat.samples(event_dat.ttl_idx2),csc_dat.dataSamples,params);
    
    ix = find(-dt/2:1/csc_dat.Fs:dt/2<=0.0004);
    %event_dat.ttlavg1(ix,:) = 0;
    %event_dat.ttlavg2(ix,:) = 0;
%%
x = sum(event_dat.ttlavg1.^2,1);
z1 = (x-mean(x))./std(x);
x = sum(event_dat.ttlavg2.^2,1);
z2 = (x-mean(x))./std(x);

del_idx = unique(sort([find(z1>2.5) find(z2>2.5)]));

event_dat.ttl_idx(del_idx) = [];
event_dat.ttl_idx1(del_idx) = [];
event_dat.ttl_idx2(del_idx) = [];

event_dat.ttlavg1(:,del_idx) = [];
event_dat.ttlavg2(:,del_idx) = [];
LogDat.RT( del_idx ) = [];
LogDat.stimID( del_idx ) = [];
   

    %%
    [LogDat.RT,LogDat.s_idx] = sort(LogDat.RT);
    [LogDat.raRT] = compute_runningAVG(LogDat.RT,6);
    LogDat.stimID = LogDat.stimID(LogDat.s_idx);

%%

    tx =-dt/2:1/csc_dat.Fs:dt/2;
    n = size(event_dat.ttlavg1,2);
    
    figure;
    subplot(2,2,1);
    a1 = gca;
    hold on;
    imagesc(tx,1:n,event_dat.ttlavg1(:,LogDat.s_idx)');
    plot([0 0],[0 n],'w','LineWidth',3);
    for it = 1:n
        plot(LogDat.raRT(it)-1,it,'.','Color',[0 0 0])
    end;
    axis xy;axis tight;
    caxis([-50 50]);
    xlim([-.5 1.5]);
    title('Stimulus locked');
    
    subplot(2,2,3);
    a2 = gca;
    hold on;
    y = mean(event_dat.ttlavg1,2);
    %plot([0 0],[min(y) max(y)],'c','LineWidth',.1);
    plot(tx,y,'r-');
    axis tight;
    xlim([-.5 1.5]);

    subplot(2,2,2);
    a1 = [a1 gca];
    hold on;
    imagesc(tx,1:n,event_dat.ttlavg2(:,LogDat.s_idx)');
    plot([0 0],[0 n],'w','LineWidth',3);
    axis xy;axis tight;
    axis xy;
    caxis([-50 50]);
    xlim([-1.5 .5]);
    title('Response locked');
    subplot(2,2,4);
    a2 = [a2 gca];
    hold on;
    y = mean(event_dat.ttlavg2,2);
    %plot([0 0],[min(y) max(y)],'c','LineWidth',.1);
    plot(tx,y,'r-');
    axis tight;
    xlim([-1.5 .5]);
    
    for it = 1:length(a1)
        xlabel(a1(it),'Time (s)');
        ylabel(a1(it),'Trial #');
    end;
    
    yl = zeros(2,2);
    for it = 1:length(a2)
        xlabel(a2(it),'Time (s)');
        ylabel(a2(it),'Amplitude (\muV)');
        yl(it,:) = get(a2(it),'YLim');
    end;
    set(a2,'YLim',[min(min(yl)) max(max(yl))]);
    %%
%     for it = 1:length(event_dat.ttl_idx1)
%         if sign(event_dat.samples(event_dat.ttl_idx1(it))-dt/2*csc_dat.Fs)==1
%             csc_dat.filtdat(floor(event_dat.samples(event_dat.ttl_idx1(it))-dt/2*csc_dat.Fs:event_dat.samples(event_dat.ttl_idx1(it))+dt/2*csc_dat.Fs)) = csc_dat.filtdat(floor(event_dat.samples(event_dat.ttl_idx1(it))-dt/2*csc_dat.Fs:event_dat.samples(event_dat.ttl_idx1(it))+dt/2*csc_dat.Fs))-mean(event_dat.ttlavg1,2);
%         end;
%     end;
%     
%     for it = 1:length(event_dat.ttl_idx2)
%         if sign(event_dat.samples(event_dat.ttl_idx2(it))-dt/2*csc_dat.Fs)==1 && (floor(event_dat.samples(event_dat.ttl_idx2(it))+event_dat.ttl_idx2(it)+dt/2*csc_dat.Fs) < length(csc_dat.dataSamples))
%             csc_dat.filtdat(floor(event_dat.samples(event_dat.ttl_idx2(it))-dt/2*csc_dat.Fs:event_dat.samples(event_dat.ttl_idx2(it))+dt/2*csc_dat.Fs)) = csc_dat.filtdat(round(event_dat.samples(event_dat.ttl_idx2(it))-dt/2*csc_dat.Fs:event_dat.samples(event_dat.ttl_idx2(it))+dt/2*csc_dat.Fs))-mean(event_dat.ttlavg2,2);
%         end;
%     end;
    %%    
    %[datafilt] = locdetrend(datafilt,Fs,[.1 .05]);
    
    %[cleanSignal noise] = CleanLineNoise(datafilt','Fs',Fs,'noiseFreq',50,'windowSize',1);
    %clear datafilt;
    %clear noise;
    %%
    event_dat.samples = event_dat.samples(event_dat.ttl_idx);
    %%
    Wp = [600 8000]./csc_dat.Fs/2;
    [b,a] = butter(4,Wp);
    
    Hd{1} = b;
    Hd{2} = a;
    
    params = [];
    params.nrNoiseTraces = 0; 
    params.detectionMethod= 1;% 1 -> from power signal, 2 threshold positive, 3 threshold negative, 4 threshold abs, 5 wavelet
    %params.detectionParams: depends on detectionMethod. 
    if params.detectionMethod==1; params.kernelSize = 18;end;
    if params.detectionMethod==4; params.scaleRanges = [];params.waveletName = '';end;
    params.alignMethod = 1;% 1-> find peak, 2->none, 3->peak of power signal, 4->peak of MTEO signal.
    params.extractionThreshold = 5;
    params.prewhiten = 0;
    params.samplingFreq = csc_dat.Fs;
    params.limit = 1000;
    
    [rawMean, filteredSignal, rawTraceSpikes,spikeWaveforms, spikeTimestamps, runStd2, upperlim, noiseTraces] = extractSpikes(csc_dat.dataSamples, Hd, params );
    %%
    params = [];
    params.Fs = csc_dat.Fs;
    dt = 0.3;
        
    [stAVG] = eta(dt,spikeTimestamps,csc_dat.dataSamples,params);
    
    figure;
    hold on;
    plot(-dt/2:1/csc_dat.Fs:dt/2,mean(stAVG,2));
%     %%
%     [AC]= compute_spiketrain_autocorr(csc_dat.tsi(spikeTimestamps)./1e6);
%     figure;
%     plot(AC.time{1},AC.avg);
    %%
    x = cell(length(event_dat.samples),1);
    spk_dat = struct;
    for it = 1:length(event_dat.samples)
        tlim = [event_dat.samples(it)-csc_dat.Fs*pre event_dat.samples(it)+csc_dat.Fs*post];
        x{it} = find(spikeTimestamps >= tlim(1) & spikeTimestamps <=tlim(2));
        spk_dat.times{it} = [csc_dat.tsi(spikeTimestamps(x{it}))./1e6]';
    end;
    
    x = x(LogDat.s_idx);
    spk_dat.times = spk_dat.times(LogDat.s_idx);
    %%
    ev = event_dat.events(event_dat.ttl_idx(LogDat.s_idx),1)./1e6;
    st = csc_dat.tsi(spikeTimestamps)./1e6;
    dt = [pre post];
    bw =0.01;
    
    [spk_dat.pSTH] = compute_pSTH(ev,st,dt,bw);    
    %%
    spk_dat.x = -dt(1):bw:post;
    spk_dat.fr = [];
    for it = 1:length(LogDat.ID)
        
        ix = find(LogDat.stimID == LogDat.ID(it));        
        
        dum = spk_dat.pSTH(ix,:);
        [spk_dat.fr(it,:)] = compute_firing_rate(dum,bw);
%         for jt = 1:length(ix)
%             ix2 = find(tx >= LogDat.RT(ix(jt)));
%             dum(jt,ix2) = 0;
%         end;
%         
%         [fr2(it,:)] = compute_firing_rate(dum,bw);
        
    end;
    %%       
    figure;
    subplot(2,1,1);
    hold on;
    plot(zeros(length(x),1),[1:length(x)],'ro');
    for it = 1:length(x)
        plot(LogDat.raRT(it)-1,it,'mo');
        plot(csc_dat.tsi(spikeTimestamps(x{it}))/1e6 - event_dat.events(event_dat.ttl_idx(LogDat.s_idx(it)),1)/1e6,it*ones(1,length(x{it})),'k.');
    end;
    axis tight;
    xlim([-pre post]);
    xlabel('Time (s)');
    ylabel('Trial #');
    
    subplot(2,1,2);
    hold on;
    y1 = mean(spk_dat.fr,1);
    plot(spk_dat.x,y1);
    axis tight;
    xlim([-pre post]);
    %%
    k1 = ceil(length(LogDat.ID)/10)+1;
    k2 = 10;
    figure;
    subplot(k1,k2,1);
    plot([0 0],[0 6],'r');
    xlabel('Time (s)');
    ylabel('Trial #');
    xlim([-pre post]);
    ylim([1 6]);
    
    for it = k2+1:length(LogDat.ID)
        subplot(k1,k2,it);
        hold on;     
        six = find(LogDat.stimID == LogDat.ID(it));
        plot([0 0],[0 length(six)],'r');
        for jt = 1:length(six)
            x = spk_dat.times{six(jt)}-event_dat.events(event_dat.ttl_idx(LogDat.s_idx(six(jt))),1)/1e6;
            for kt = 1:length(x)
                plot([x(kt) x(kt)],[jt-1 jt],'k-')
            end;
        end;
        axis tight;
        xlim([-pre post]);
        str = ['pic#:',num2str(LogDat.ID(it))];
        title(str);
    end;
%     %%
%     k1 = round(length(fit)/10);
%     k2 = 10;
%     figure;
%     for it = 1:length(fit)
%         subplot(k1,k2,it);
%         lfplot(fit{it});
%         lfband(fit{it});
%         xlim([-pre post]);
%     end;
%     %%
%     figure;
%     subplot(5,1,1:4);imagesc(-pre:1/32000:post,1:size(datafilt,2),datafilt');
%     subplot(515);plot(-pre:1/32000:post,mean(datafilt,2),'b-');
end;