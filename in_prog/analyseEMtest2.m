function analyseEMtest2(pID,expMode,sesh,nWorkers)

if nargin ==0
    pID                     = 'P07';
    expMode                 = 'fVSpEM';
    sesh                    = '2017-04-30_17-35-52';
    nWorkers                = 12;
end;

%%
addpath('/home/rouxf/tbx/releaseDec2015/');
addpath(genpath('/home/rouxf/tbx/osort-v3-rel/'));
addpath(genpath('/home/rouxf/tbx/wave_clus/'));
addpath(genpath('/home/rouxf/tbx/chronux_2_11/'));

addpath('/home/rouxf/prj/Bham/code/mcode/params/');
addpath(genpath('/home/rouxf/prj/Bham/code/mcode/helper/'));
addpath('/home/rouxf/prj/Bham/code/mcode/utils/spikes/');

%%
Fs = 32e3;
[par] = set_parameters_Bham(Fs);

%%
[rpath] = ['/media/rouxf/rds-share/Archive/MICRO/',pID,'/',expMode,'/'];

x = dir(rpath);
x(1:2) = [];

[p2d] = [rpath,sesh,'/'];
[CSCfiles] = dir([p2d,'*.ncs']);

[savepath] = ['/home/rouxf/out/',expMode,'/',sesh,'/'];

%[savepath] = p2d;
%savepath = [savepath(1:regexp(savepath,'Archive')-1),'iEEG_DATA',savepath(regexp(savepath,'Archive/')+7:end)];

chck = dir(savepath);
if isempty(chck)
    mkdir(savepath);
end;

chck = dir([savepath,'lfp_dat/']);
if isempty(chck)
    mkdir([savepath,'lfp_dat/']);
end;

chck = dir([savepath,'spike_dat/']);
if isempty(chck)
    mkdir([savepath,'spike_dat/']);
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

if ~isempty(EVfile)
    [TimeStamps, ttls, Hdr] = Nlx2MatEV_v3( [p2d,EVfile.name], FieldSelection, ExtractHeader, ExtractMode, ModeArray );
    
    [events] = zeros(size(ttls,2),2);
    events(:,1) = TimeStamps';
    events(:,2) = ttls';
end;

%%
[savename] = [pID,'_',sesh,'_LOGdat.mat'];

[readme] = {'TimeStamps','events','ttls','Hdr'};

save([savepath,savename],'TimeStamps','events','ttls','Hdr','readme');
    
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
    parpool(nWorkers,'SpmdEnabled',false);
end;

%%
parfor it = 1:2%1:length(CSCfiles)
    
    fprintf([num2str(it),'/',num2str(length( CSCfiles) )]);
    
    %%
    dum = CSCfiles(it).name;
    dum(regexp(dum,'_')) = [];
    dum(regexp(dum,'CSC'):regexp(dum,'CSC')+2) = [];
    dum(regexp(dum,'.ncs'):end) = [];
    chanLab = dum;
    
    %%
    
    [timestamps, dataSamples,hdr] = Nlx2MatCSC_v3([p2d,CSCfiles(it).name], FieldSelection, ExtractHeader, ExtractMode, []);
    [timestamps] = timeStampinter(timestamps);
    
    if length(dataSamples(:)) ~= length( timestamps )
        error('timestamp interpolation error detected!');
    end;
    
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
    [b,a] = butter(4,300./(Fs/2),'low');% apply low-pass for LFP
    [LFPsig] = filtfilt(b,a,dataSamples);
    
    %%
    [LFPsig,~] = CleanLineNoise( LFPsig ,'Fs', Fs , 'noiseFreq', 50,'windowSize',1);        
    
    %%
    [lfpTime] = 0:1/Fs:(size(LFPsig,2)-1)/Fs;
    
    if length(lfpTime) ~= length( timestamps )
        error('timestamp interpolation error detected!');
    end;
    
    %%
    [savename] = [pID,'_',sesh,'_',chanLab,'_LFPdat.mat'];
    
    readme = {'timeStamps','lfpDat','Fs','lfpTime','chanLab'};
    par_save([[savepath,'lfp_dat/'],savename],[],{timestamps,LFPsig,Fs,lfpTime,chanLab,readme});

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
    [savename] = [pID,'_',sesh,'_',chanLab,'_SPKdat.mat'];
    
    readme = {'sortedSpikes','wltCoeffs','timestamps'};
    par_save([[savepath,'spike_dat/'],savename],[],{sortedSpikes,wltCoeffs,timestamps,readme});
    
    %%        
    fprintf('\n');
    
end;

%%
delte(gcp);
exit;