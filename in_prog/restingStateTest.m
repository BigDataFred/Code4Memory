%% set Matlab path environment
addpath('/home/rouxf/tbx/releaseDec2015/');
addpath(genpath('/home/rouxf/tbx/osort-v3-rel/'));
addpath(genpath('/home/rouxf/tbx/wave_clus/'));
addpath(genpath('/home/rouxf/tbx/chronux_2_11/'));
addpath('/home/rouxf/tbx/fieldtrip-20170618/');
ft_defaults;

addpath('/home/rouxf/prj/Bham/code/mcode/params/');
addpath('/home/rouxf/prj/Bham/code/mcode/helper/');
addpath('/home/rouxf/prj/Bham/code/mcode/utils/spikes/');
addpath('/home/rouxf/prj/Bham/code/mcode/utils/LFP/');
addpath('/home/rouxf/prj/Bham/code/mcode/project_EM/');

%% recruit workers
% Create a parallel pool if none exists
if isempty(gcp('nocreate'))
    parpool(36,'SpmdEnabled',false);
end
    
%% set data labels
pID = 'P07';
expMode = 'fVSpEM';

%% set data parameters
Fs =32e3;
timeStampsPerSec = 1/Fs*1e6;

[par] = set_parameters_Bham(Fs);

%% parameters to read the TTL data    
FieldSelection1(1) = 1;%timestamps
FieldSelection1(2) = 0;% EventIDs
FieldSelection1(3) = 1;%TTLs
FieldSelection1(4) = 0;% Extras
FieldSelection1(5) = 0;% Event strings
    
%% parameters to read the CSC data
FieldSelection2(1) = 1;%timestamps
FieldSelection2(2) = 0;
FieldSelection2(3) = 0;%sample freq
FieldSelection2(4) = 0;
FieldSelection2(5) = 1;%samples
    
%% search for existing data
%[rpath] = '/media/rouxf/rds-share/Archive/MICRO/P09/Spontaneous/';

%[rpath] = '/media/rouxf/rds-share/Archive/MICRO/P07/cnEM/';

[rpath] = '/home/rouxf/in/';
[savepath] = ['/home/rouxf/out/',expMode,'/'];

sesh = dir(rpath);
sesh(1:2) = [];

[chck] = regexp({sesh(:).name},'\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}');

sel = [];
for it = 1:length(chck)
    if ~isempty(chck{it})
        sel = [sel it];
    end;
end;
sesh = sesh(sel);

%% loop over recording sessions
for seshSel = 1:length(sesh)
    
    [p2d] = [rpath,sesh(seshSel).name,'/'];
    %[CSCfiles] = dir([p2d,'CSC_postHippR*.ncs']);
    [CSCfiles] = dir([p2d,'*.ncs']);
    
    %% search the Logfile
    p2logDat = ['/media/rouxf/rds-share/iEEG_Data/MICRO/',pID,'/',expMode,'/',sesh(seshSel).name,'/log_dat/'];
    logFile = dir([p2logDat,'*_Log*ile_EMtask_LogDat.mat']);
        
    %% load the logfile data    
    [logDat] = load([p2logDat,logFile.name]);    
    
    %% read the TTL event data    
    [EVfile] = dir([p2d,'*.nev']);
    
    [TimeStamps, ttls, Hdr] = Nlx2MatEV_v3( [p2d,EVfile.name], FieldSelection1, 1, 1, [] );
    
    [events] = zeros(size(ttls,2),2);
    events(:,1) = TimeStamps';
    events(:,2) = ttls';
    firstTstamp = events(1,1);
    
    events(:,1) = (events(:,1)-firstTstamp)./1e6;
        
    selIx = find(events(:,2)==7);
    events = events(selIx,:);
   
    %%
    sortedSpikes = cell(1,length(CSCfiles));
    wltCoeffs = cell(1,length(CSCfiles));
    wvf = cell(1,length(CSCfiles));
    chanLab = cell(1,length(CSCfiles));
    LFPsig = cell(1,length(CSCfiles));
    
    %%
    %for it = 5%1:length(CSCfiles)
    parfor it = 1:length(CSCfiles)
        
        tic;
        fprintf([num2str(it),'/',num2str(length( CSCfiles) )]);
        
        %%
        dum = CSCfiles(it).name;
        dum(regexp(dum,'_')) = [];
        dum(regexp(dum,'CSC'):regexp(dum,'CSC')+2) = [];
        dum(regexp(dum,'.ncs'):end) = [];
        chanLab{it} = dum;
        
        %%
        [timestamps, dataSamples,hdr] = Nlx2MatCSC_v3([p2d,CSCfiles(it).name], FieldSelection2, 1, 1, []);
        
        chck = regexp(hdr,'ADBitVolts');
        selIdx = [];
        for jt = 1:length(chck);selIdx(jt) = ~isempty(chck{jt});end;
        selIdx = find(selIdx~=0);
        scalef = str2double(hdr{selIdx}(min(regexp(hdr{selIdx},'\d')):end));
        
        %     chck = regexp(hdr,'SamplingFrequency');
        %     selIdx = [];
        %     for jt = 1:length(chck);selIdx(jt) = ~isempty(chck{jt});end;
        %     selIdx = find(selIdx~=0);
        %     Fs = str2double(hdr{selIdx}(min(regexp(hdr{selIdx},'\d')):end));
        
        %% flatten
        [dataSamples] = double(dataSamples(:))';
        [dataSamples] = dataSamples.*scalef.*1e6;
        
        %%
        LFPsig{it} = dataSamples;
        
        %% spike detection
        par2 = par;
        par2.fnamespc = [par2.fnamespc,num2str(it)];
        par2.fname_in = [par2.fname_in,num2str(it)];
        
        [~,spikeWaveforms,~,spikeTimestamps,~,~,noise_std_undetect] = amp_detect(dataSamples,par2);
        
        %%
        waveclus.spikes                     = spikeWaveforms;
        waveclus.index                      = spikeTimestamps;
        
        [dim] = size(spikeWaveforms);
        
        %% do spike sorting
        par2.filename = [CSCfiles(it).name];
        
        if dim(1) > dim(2)% the number of spike events must larger than the number of samples in the waveform
            [dum,wltCoeffs{it}] = doSpikeSorting_waveclus( waveclus , par2 );
        end;
        dum.sd                = noise_std_undetect;
        
        %% force cluster-membership
        selIx1 = find(dum.assignedCluster ==0);
        selIx2 = find(dum.assignedCluster ~=0);
        [class_out] = force_membership_wc(wltCoeffs{it}(selIx2,:), dum.assignedCluster(selIx2), wltCoeffs{it}(selIx1,:), par);
        dum.assignedCluster(selIx1(class_out~=0)) = class_out(class_out~=0);
        
        %% delete noise clusters
        delIx = find(dum.assignedCluster==0);
        dum.assignedCluster( delIx ) = [];
        dum.newSpikeTimes( delIx ) = [];
        dum.wavf(delIx,:) = [];
        
        %%
        sortedSpikes{it} = dum;
        
        %%
        waveclus = [];
        
        %% interpolate LFP
        [dum] = interpLFP(LFPsig{it},spikeTimestamps,[0.002 0.006],Fs,'linear');
        
        %% lowpass filter LFP
        [b,a] = butter(4,[300]./(Fs/2),'low');% apply low-pass for LFP
        [dum] = filtfilt(b,a,dum);
        
        %%
        [dum,noise] = CleanLineNoise(dum,'Fs',Fs,'noiseFreq',50,'windowSize',1);
        
        %     %% notch filter LFP
        %     [b,a] = butter(2,[49.5 50.5]./(Fs/2),'stop');% apply band-stop for LFP
        %     dum2 = filtfilt(b,a,dum);
        %
        %     [b,a] = butter(2,[99.5 100.5]./(Fs/2),'stop');% apply band-stop for LFP
        %     dum = filtfilt(b,a,dum);
        
        %%
        [LFPsig{it}] = dum;
        dum = [];
        
        %%
        fprintf('\n');
        toc;
        
    end;
    fprintf('Finished raw data import & preprocessing, moving on\n');
    
    %% time vector for LFP
    [lfpTime] = [1:size(LFPsig{1},2)]./Fs;
    
    % savename = [pID,'_',expMode,'_',sesh(seshSel).name,'_lfpDataStimLocked.mat'];
    % save([savepath,savename],'LFPsig','lfpTime','Fs','chanLab','-v7.3');
    %
    % savename = [pID,'_',expMode,'_',sesh(seshSel).name,'_spkDataStimLocked.mat'];
    % save([savepath,savename],'sortedSpikes','chanLab','-v7.3');
    
    %%    
    [RTs]    = str2double(logDat.LogDat1.log(:,end));        
    [trlENC] = extract_Enc_trl(logDat.LogDat1);
    
    %% indexes of trials that were later remebered and forgotten
    [hitIdx] =logDat.ix{4};%HITS
    [missIdx] = [logDat.ix{5};logDat.ix{6}];%MISSES
    
    %% delete trials with RTs < 4 sec
    delIx = find(RTs <=4);
    
    trlENC(delIx) = [];        
    RTs(delIx) = [];
    hitIdx(find(ismember(hitIdx,delIx))) = [];
    missIdx(find(ismember(missIdx,delIx))) = [];
    
    %% readout events corresponding to Encoding
    [events] = events(trlENC,:);
    
    %% compute the trl matrix
    pre = 5;
    post = 7;
    [trl] = zeros(size(events,1),3);
    for it = 1:size(events,1)
        ix = find( lfpTime >= (events(it,1)-pre) & lfpTime < (events(it,1)+post) );
        [trl(it,:)] = [min(ix) max(ix) -pre*Fs];
    end;
    
    %% segment the LFP data
    ntrl = size(trl,1);
    nsmp = length(trl(1,1):trl(1,2));
    trlTime = -pre:1/Fs:post-(1/Fs);
    
    LFPseg = cell( 1 , length(LFPsig) );
    for jt = 1:length(LFPsig)
        
        fprintf([num2str(jt),'/',num2str(length(LFPsig))]);
        
        x = zeros(size(trl,1),(pre+post)*Fs);
        for it = 1:size(trl,1)
            x(it,:) = LFPsig{jt}(trl(it,1):trl(it,2));
        end;
        LFPseg{jt} =  x';
        
        fprintf('\n');
    end;
    clear LFPsig trl;
    
    %%
    % savename = [pID,'_',expMode,'_',sesh(seshSel).name,'_lfpDataStimLockedSegmented.mat'];
    % save([savepath,savename],'LFPseg','trlTime','Fs','chanLab','-v7.3');
    
    %% segment the spike data
    parfor jt = 1:length(sortedSpikes)
        lfpTime;
        [spkTime] = lfpTime(sortedSpikes{jt}.newSpikeTimes);
        sortedSpikes{jt}.SpikeTimesSeg = [];
        sortedSpikes{jt}.trl = [];
        sortedSpikes{jt}.assignedClusterSeg = [];
        for it = 1:size(events,1)
            
            [ix] = find( spkTime >= (events(it,1)-pre) & spkTime <= (events(it,1)+post) );
            
            sortedSpikes{jt}.SpikeTimesSeg= [sortedSpikes{jt}.SpikeTimesSeg spkTime(ix)-events(it,1)];
            sortedSpikes{jt}.trl = [sortedSpikes{jt}.trl it*ones(1,length(ix))];
            sortedSpikes{jt}.assignedClusterSeg = [sortedSpikes{jt}.assignedClusterSeg sortedSpikes{jt}.assignedCluster(ix)];
        end;
    end;
    
    %%
    savename = [pID,'_',expMode,'_',sesh(seshSel).name,'_spkDataStimLockedSegmented.mat'];
    save([savepath,savename],'sortedSpikes','Fs','chanLab','RTs','trlENC','hitIdx','missIdx','-v7.3');
    clear sortedSpikes;
    
    %% use fieldtrip to downsample LFP data
    parfor it = 1:length( LFPseg )
        
        dum = [];
        dum.label = {'dumChan'};
        dum.trial = cell(1,size(LFPseg{it},2));
        dum.time = cell(1,size(LFPseg{it},2));
        for jt = 1:size( LFPseg{it},2)
            dum.trial{jt} = LFPseg{it}(:,jt)';
            dum.time{jt} = -pre:1/Fs:post-1/Fs;
        end;
        
        cfg                     = [];
        cfg.resamplefs          = 1e3;
        
        [dum] = ft_resampledata( cfg, dum );
        
        cfg                     = [];
        cfg.keeptrials          = 'yes';
        
        [dum] = ft_timelockanalysis( cfg , dum );
        
        LFPseg{it} = [];
        LFPseg{it} = squeeze(dum.trial)';
        
    end;
    Fs = 1e3;
    trlTime = -pre:1/Fs:post-(1/Fs);
    
    % %%
    % parfor it = 1:length( LFPseg )
    %     for jt = 1:size(LFPseg{it},2)
    %
    %         [dum,noise] = CleanLineNoise(LFPseg{it}(:,jt)','Fs',1e3,'noiseFreq',50,'windowSize',1);
    %         LFPseg{it}(:,jt) = dum';
    %     end;
    % end;
    
    %%
    savename = [pID,'_',expMode,'_',sesh(seshSel).name,'_lfpDataStimLockedSegmenteddownsampled.mat'];
    save([savepath,savename],'LFPseg','trlTime','Fs','chanLab','RTs','trlENC','hitIdx','missIdx','-v7.3');
    
    %% compute average LFP across BF electrodes
    chck = regexp(chanLab,'\d{1,8}');
    chck  = [chck{:}];
    
    BF = cell(1,length(chanLab));
    parfor it = 1:length(chanLab)
        BF(it) = {chanLab{it}(1:chck(it)-1)};
    end;
    BFid = unique(BF);
    
    [LFPavg] = cell( 1 , length(BFid) );
    parfor it = 1:length(BFid)
        LFPavg{it} = zeros(size(LFPseg{1}));
        
        ix = find(strcmp(BF,BFid(it)));
        for jt = 1:length(ix)
            LFPavg{it} = LFPavg{it} + LFPseg{ix(jt)};
        end;
        LFPavg{it} = LFPavg{it}./length(ix);
    end;
    clear LFPseg;
    chanLab = BFid;
       
    %%
    savename = [pID,'_',expMode,'_',sesh(seshSel).name,'_lfpDataStimLockedSegmentedAVGdownsampled.mat'];
    save([savepath,savename],'LFPavg','trlTime','Fs','chanLab','RTs','trlENC','hitIdx','missIdx','-v7.3');
    
end;
%% shut down parpool
delete(gcp);

%%
movingwin1 = [1 0.01];
T = movingwin1(1);
W = 2;
TW = T*W;

params1                  =[];
params1.pad              = 0;
params1.Fs               = 1e3;
params1.fpass            = [0 30];
params1.tapers           = [TW 2*TW-1];
params1.trialave         = 1;

movingwin2 = [0.5 0.01];
T = movingwin2(1);
W = 10;
TW = T*W;

params2                  =[];
params2.pad              = 0;
params2.Fs               = 1e3;
params2.fpass            = [30 150];
params2.tapers           = [TW 2*TW-1];
params2.trialave         = 1;

Sl = cell(1,length(LFPavg));
Sh = cell(1,length(LFPavg));
for it = 1:length(LFPavg)
    fprintf([num2str(it),'/',num2str(length(LFPavg))]);
    [Sl{it},tl,fl] = mtspecgramc( gradient(LFPavg{it}')', movingwin1, params1 );
    [Sh{it},th,fh] = mtspecgramc( gradient(LFPavg{it}')', movingwin2, params2 );
    fprintf('\n');
end;

%%
for jt = 1:length(sortedSpikes)
    [spkTime] = lfpTime(sortedSpikes{jt}.newSpikeTimes);
    sortedSpikes{jt}.SpikeTimesSeg = [];
    sortedSpikes{jt}.trl = [];
    sortedSpikes{jt}.assignedClusterSeg = [];
    for it = 1:size(events,1)
        
        [ix] = find( spkTime >= (events(it,1)-pre) & spkTime <= (events(it,1)+post) );
        
        sortedSpikes{jt}.SpikeTimesSeg= [sortedSpikes{jt}.SpikeTimesSeg spkTime(ix)-events(it,1)];
        sortedSpikes{jt}.trl = [sortedSpikes{jt}.trl it*ones(1,length(ix))];
        sortedSpikes{jt}.assignedClusterSeg = [sortedSpikes{jt}.assignedClusterSeg sortedSpikes{jt}.assignedCluster(ix)];
    end;
end;



%%
c = 0;
ts = {};
for it = 1:length( sortedSpikes )
    
    cID = unique(sortedSpikes{it}.assignedClusterSeg);
    for jt = 1:length(cID)
        c = c+1;
        [ix] = find(sortedSpikes{it}.assignedCluster == cID(jt));
        [ts{c}] = lfpTime(sortedSpikes{it}.newSpikeTimes).*1e3;
        
    end;
end;

%%
[dt] = lfpTime.*1e3;
c = 0;
xc = cell(1,length(ts)-1);
ix =  setdiff(26:40,34);
for it = 34%:length( ts )
    
    fprintf([num2str(it),'/',num2str(length( ts ))]);
    %ix = setdiff(1:length( ts ), it);
    x1 = ts{it};
    [x1,~] = hist(x1,dt);
    
    for jt = 1:length( ix )
        c = c+1;
        x2 = ts{ix(jt)};
        [x2,~] = hist(x2,dt);
        
        xc{c} = xcorr(x1,x2,500);
    end;
    fprintf('\n');
    
end;

%%
T = 12;
W = 1;
TW = T*W;
params1                     = [];
params1.Fs                  = 1e3;
params1.pad                 = 0;
params1.tapers              = [TW 2*TW-1];
params1.fpass               = [0 30];
params1.trialave            = 1;

T = 12;
W = 2;
TW = T*W;
params2                     = [];
params2.Fs                  = 1e3;
params2.pad                 = 0;
params2.tapers              = [TW 2*TW-1];
params2.fpass               = [0 30];
params2.trialave            = 1;

c = 0;
isiH = [];
isiP = [];
wvf  = {};
xc = [];
S = {};
C = {};
R = {};
chIx = [];
for it = 1:length(sortedSpikes)
    
    fprintf([num2str(it),'/',num2str(length(sortedSpikes))]);
    
    cID = unique(sortedSpikes{it}.assignedCluster);
    for jt = 1:length(cID)
        
        if cID(jt) ~=0
            
            c = c+1;
            
            chIx(c,:) = [it cID(jt)];
            
            [ix] = find(sortedSpikes{it}.assignedCluster == cID(jt));
            
            wvf{c} = sortedSpikes{it}.wavf(ix,:);
            ts = lfpTime(sortedSpikes{it}.newSpikeTimes(ix)).*1e3;
            
            ix2 = find(sortedSpikes{it}.assignedClusterSeg == cID(jt));
            trl = sortedSpikes{it}.trl(ix2);
            ts2 = sortedSpikes{it}.SpikeTimesSeg(ix2).*1e3;
            
            ts = ts-ts(1);
            dx = diff(ts);
            isiP(c) = (length(find(dx<=3))/length(dx))*100;
            
            dt = 0:1:500;
            dx(dx>500) = [];
            [isiH(c,:),~] = hist(dx,dt);
            
            dt = min(ts):1:max(ts);%lfpTime.*1e3;
            dum = ts;
            dum(dum<dt(1) & dum >dt(end)) = [];
            
            [x,~] = hist(dum,dt);
            xc(c,:) = xcorr(x,1.5e3);
            
            dt = -pre*1e3:1:(post)*1e3-1;
            data2 = zeros(size(LFPavg{1}));
            data1 = cell(1,length(LFPavg));
            for kt = 1:size(LFPavg{1},2)
                
                ix3 = find(trl == kt);
                ts3 = ts2(ix3);
                
                ts3(ts3<dt(1)) = [];
                ts3(ts3>dt(end)) = [];
                
                [n,~] = hist(ts3,dt);
                data2(:,kt) = n;
                
                spkix = find(n==1);
                for lt = 1:length(LFPavg)
                    if ~isempty(spkix)
                        [data1{lt}(:,kt)] = interpLFP(LFPavg{lt}(:,kt),spkix,[0.002 0.006],Fs,'linear');
                    else
                        [data1{lt}(:,kt)] = LFPavg{lt}(:,kt);
                    end;
                end;
                
            end;
            
            [S{c},f1,R{c}] = mtspectrumpb(data2,params1);
            
            for lt = 1:length(data1)
                [C{c,lt}] = coherencycpb(data1{lt},data2,params2);
            end;
            
        end;
        
    end;
    fprintf('\n');
    
end;

%%
Clr = [0 0 0;...
    173 255 47;...
    0 255 255;...
    255 165 0;...
    255 69 0;...
    139 0 0];

Clr = Clr./255;

for it = 26:40%1:size(isiH,1)
    
    figure;
    subplot(231);
    hold on;
    plot(linspace(0,2,64),wvf{it},'r');
    plot(linspace(0,2,64),mean(wvf{it},1),'k');
    title([chanLab{chIx(it,1)},' c#:', num2str(chIx(it,2))]);
    xlabel('Time [ms]');
    ylabel('[\muV]');
    set(gca,'XTick',[0 2]);
    
    subplot(232);
    bar(0:500,isiH(it,:));
    axis tight;
    xlim([-1 501]);
    title([num2str(round(isiP(it)*10)/10),'% < 3ms']);
    xlabel('ISI [ms]');
    ylabel('Count');
    set(gca,'XTick',[0 250 500]);
    
    subplot(234);
    plot(1:1500,xc(it,1502:end));
    axis tight;
    xlabel('Lag [ms]');
    ylabel('Coincidences');
    set(gca,'XTick',[0 0.5 1].*1e3);
    
    subplot(235);
    plot( f1 , S{it} );
    axis tight;
    title([num2str(round(R{it}*10)/10),'Hz']);
    xlabel('Freq. [Hz]');
    ylabel('Spike-Power [a.u.]');
    set(gca,'XTick',[5 15 25]);
    
    subplot(236);
    hold on;
    sIx = [find(strcmp(BF{chIx(it,1)},BFid)) setdiff(1:length(BFid),find(strcmp(BF{chIx(it,1)},BFid)))];
    for jt = 1:size(C,2)
        plot( f1 , C{it,jt} ,'Color', Clr(sIx(jt),:));
    end;
    axis tight;
    xlabel('Freq. [Hz]');
    ylabel('Spike-LFP Coherence [a.u.]');
    set(gca,'XTick',[5 15 25]);
    
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

for it = 1:length(sortedSpikes)
    
    fprintf([num2str(it),'/',num2str(length( CSCfiles) )]);
    
    [cluID] = unique(sortedSpikes{it}.assignedCluster);
    
    for ct = 1:length(cluID)
        
        if ~isempty(cluID) && cluID(ct) ~=0
            if ~isempty(sortedSpikes{it}.newSpikeTimes)
                
                c = c+1;
                
                cluInf(c,1) = it;
                cluInf(c,2) = cluID(ct);
                
                ix = find( sortedSpikes{it}.assignedCluster == cluID(ct) );
                
                ts{c} = lfpTime( sortedSpikes{it}.newSpikeTimes( ix ) );
                ts{c} = ts{c}.*1e3;
                
                wvf{c} = sortedSpikes{it}.wavf(ix,:);
                
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
selIdx = find(ismember(spikeTimestamps,idx));
wvf = spikeWaveforms(selIdx,:);
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