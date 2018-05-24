%%
clear;clc;

%% open parpool
%if isempty(gcp('nocreate'))
%    %parpool(length(params.CSC_idx),'SpmdEnabled',false);
%    parpool(4);
%end;

%% set the path def
restoredefaultpath;
addpath(genpath('/media/rouxf/rds-share/Fred/code/mcode/custom/'));% custom code
addpath(genpath('/media/rouxf/rds-share/Common/releaseDec2015/'));%needed to read Nlx data
addpath(genpath('/media/rouxf/rds-share/Common/osort-v3-rel/'));%needed for spike detection & sorting
addpath(genpath('/media/rouxf/rds-share/Common/chronux_2_11/'));% needed for spectral analysis
addpath(genpath('/media/rouxf/rds-share/Common/wave_clus/'));%needed for spike detection & sorting
fn = dir('/media/rouxf/rds-share/Common/fieldtrip-*');
addpath(['/media/rouxf/rds-share/Common/',fn.name]);% spectral analysis, signal processing, spike detection
ft_defaults;

%% set the data-paths directions
params = [];
params.CSC_idx = [ 25:36 ];
params.pID = 'P08';
params.sesh = [ 2 ];
params.p2Nlxdata = ['/media/rouxf/rds-share/Archive/MICRO/',params.pID,'/Tunings/'];%path 2 Nlx data
params.savepath = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',params.pID,'/Tunings/'];%path for saving
params.rpath = '/home/rouxf/';

params.spike_tresh = 4;%threshold for spike detection

params.mode         = 'stimlocked';%'resplocked';
params.pre          = 2;% baseline time range
params.post         = 2.5;% params.post stim time range

% set flags
params.makeLFP      = 0;
params.makeSpikes   = 1;
params.makeMUA      = 0;

chck = dir( params.savepath );
if isempty( chck )
    mkdir( params.savepath )
end;

%%
params.ref = [];

params.ref = zeros(1,length(params.CSC_idx));

%% create the session labels of the Nlx data
params.Nlxdat = {};
[params.Nlxdat] = set_Nlx4tuning(params.p2Nlxdata,params.pID)

params.p2Nlxdata = [params.p2Nlxdata,params.Nlxdat{params.sesh},filesep];
chans = dir([params.p2Nlxdata,'*.ncs']);

params.CSClabels = {};
for it = 1:length(chans)
    params.CSClabels(it) = {chans(it).name};
end;
params.CSClabels = params.CSClabels';

chck = dir([params.savepath,params.Nlxdat{params.sesh}]);
if isempty(chck)
    mkdir([params.savepath,params.Nlxdat{params.sesh}]);
end;

%% make folders
chck = dir([params.savepath,params.Nlxdat{params.sesh}, filesep, 'log_dat']);
if isempty(chck)
    mkdir([params.savepath,params.Nlxdat{params.sesh}, filesep, 'log_dat']);
end;
    
chck = dir([params.savepath,params.Nlxdat{params.sesh}, filesep, 'spike_dat']);
if (params.makeSpikes ==1 )
    if isempty(chck)
        mkdir([params.savepath,params.Nlxdat{params.sesh}, filesep, 'spike_dat']);
    end;
end;

chck = dir([params.savepath,params.Nlxdat{params.sesh}, filesep, 'lfp_dat']);
if (params.makeLFP ==1 )
    if isempty(chck)
        mkdir([params.savepath,params.Nlxdat{params.sesh}, filesep, 'lfp_dat']);
    end;
end;

chck = dir([params.savepath,params.Nlxdat{params.sesh}, filesep, 'mua_dat']);
if (params.makeMUA ==1 )
    if isempty(chck)
        mkdir([params.savepath,params.Nlxdat{params.sesh}, filesep, 'mua_dat']);
    end;
end;

%% create the session labels of the log data
params.Tlogfile = {};
[params.Tlogfile] = set_log4tuning(['/media/rouxf/rds-share/Archive/MICRO/',params.pID,'/log/Tuning/'],params.pID)

%% this is a patch written for P02 - FIX ME
ix      = regexp(params.p2Nlxdata,'/');
NLXts     = params.p2Nlxdata(ix(end-1):end);
NLXts(regexp(NLXts,'/')) = [];
sel = params.sesh;

switch params.pID
    case 'P02'
        d       = NLXts(regexp(NLXts,'_')-2:regexp(NLXts,'_')-1);
        d(regexp(d,'-')) = [];
        d = str2double( d );
        
        h       = NLXts(regexp(NLXts,'_')+1:end-1);
        h(regexp(h,'-')) = [];
        h = h(1:4);
        h = str2double( h );
        
        d2 = []; h2 = [];
        for it = 1:length(params.Tlogfile)
            x = params.Tlogfile{it};
            
            dum = params.Tlogfile{it}(end-11:end-4);
            dum(regexp(dum,'_')) = [];
            dum = dum(1:4);
            
            d2(it) = str2double( x(min(regexp(x,'-'))-2:min(regexp(x,'-'))-1) );
            h2(it) = str2double ( dum );
        end;
        sel = find(d-d2 == 0);
        [v,ix] = min(abs(h2(sel)-h));
        sel = sel(ix+1);
end;

%% get the logfile data
f = [];
f.p2logf = ['/media/rouxf/rds-share/Archive/MICRO/',params.pID,'/log/Tuning/'];
f.logf = params.Tlogfile{sel};
p.ncols = 8;

try
    [LogDat] = getNewLogDataTune(f,p);
catch
    % this is a patch written for P05
    [LogDat] = fix_corrupted_Logfile(f,8,6,66);
end;

if length(unique(LogDat.dat(find(LogDat.stimID == LogDat.ID(1)),2))) >1
    error('wrong stimulus event assignment');
end;

[~,fn,~] = fileparts(params.Tlogfile{sel});

outname = [ fn , '_LogDat.mat' ];

save( [ params.savepath , params.Nlxdat{params.sesh},  filesep, 'log_dat' , filesep, outname ] , 'LogDat' );

%% get the event file data
[dum] = getRawTTLs([params.p2Nlxdata,'Events.nev']);
event.timestamp = dum(:,1);
event.value     = dum(:,2);
%[event] = ft_read_event([params.p2Nlxdata,'Events.nev']);

%% get onset events
ttl_idx = find([event.value(:)] == 7);
ttl_idx(find([event.value(ttl_idx+1)]>7)) = [];% remove false bits
ttl_idx = ttl_idx';

dx = diff(ttl_idx);
dx = zscore(dx);
ttl_idx(find(dx >2.5));

%del_idx = find(diff(ttl_idx)~=median(diff(ttl_idx)));
%ttl_idx(end) = [];

if LogDat.n ~= length(ttl_idx)
    error('number of events is out of range');
end;

%% find the offset events
k = 0;
ttl_idx2 = zeros(length(ttl_idx),1);
for it = 1:length(ttl_idx)
    for jt = 1:2
        if event.value(ttl_idx(it)+jt) ==0
            k = k+1;
            ttl_idx2(k) = ttl_idx(it)+jt;
            break;
        end;
    end;
end;

if length(ttl_idx2) ~= length(ttl_idx)
    error('number of events must be equal');
end;

%% extract Fs and ADV

FieldSelection(1) = 1;%timestamps
FieldSelection(2) = 1;
FieldSelection(3) = 1;%sample freq
FieldSelection(4) = 1;
FieldSelection(5) = 1;%samples
ExtractHeader = 1;

ExtractMode = 1; % 2 = extract record index range; 4 = extract timestamps range.
%ModeArray(1)=[];
%ModeArray(2)=[];

[params.timestamps, ChanN, params.sampleFreq, nrADSamples, dataSamples , headerInfo] = Nlx2MatCSC_v3([params.p2Nlxdata,params.CSClabels{1}], FieldSelection, ExtractHeader, ExtractMode, []);
dim = size(dataSamples);

%flatten
dataSamples=dataSamples(:);
nrSamples = length(dataSamples);

%%
%[params.timestamps,~,nrSamples,params.sampleFreq,~,headerInfo] = getRawCSCTimestamps( [params.p2Nlxdata,params.CSClabels{1}] );% read in header-info

params.hdr = [];
params.hdr.Fs = unique(params.sampleFreq);

[idx] = getNlxHeaderParam(headerInfo,'ADBitVolts');% conversion factor
params.hdr.scalef = str2double(headerInfo{idx}(12:end));

%[idx] = getNlxHeaderParam(headerInfo,'AcqEntName');% CSC name
%params.hdr.label = headerInfo{idx}(13:end);
%params.hdr.nChans = length(params.hdr.label);
params.hdr.label  = [];
params.hdr.nChans = [];

params.hdr.FirstTimeStamp = int64(params.timestamps(1));% first ts
params.hdr.nSamples = nrSamples;
params.hdr.TimeStampPerSample = (params.timestamps(2)-params.timestamps(1))/512;

%% convert event ts to samples
for it = 1:length([event.timestamp(:)])
    event.samples(it) = (event.timestamp(it) - double(params.hdr.FirstTimeStamp))./params.hdr.TimeStampPerSample + 1;    
end;

%%
del_idx = find([event.samples(ttl_idx)] > nrSamples);

if ~isempty( del_idx )          
    
    params.interrupted =1;
    
    datat   = (params.timestamps(:));
    ttlt    = (event.timestamp(:));
    
    dx = diff(datat);
    ix = find(dx == max(dx));
    
    f = datat(ix+1);
    
    sel = find(ttlt >= f);  
    event.timestamp = event.timestamp(sel);
    event.value = event.value(sel);
    event.samples = [];
    ttl_idx = find(event.value == 7);
    
    datat   = (params.timestamps(ix+1:end));
    ttlt = ttlt(sel);
        
    params.hdr.FirstTimeStamp = uint64(f);
    params.timestamps = datat;
    nrSamples = length(params.timestamps)*512
    
    event.timestamp = ttlt;
    for it = 1:length([event.timestamp(:)])
        event.samples(it) = (event.timestamp(it) - double(params.hdr.FirstTimeStamp))./params.hdr.TimeStampPerSample + 1;
    end;
    params.interrupted = 1;
else
    params.interrupted = 0;
end;

%% make a few sanity checks
if unique([event.value(ttl_idx)]) ~=7
    error('trigger values are out of range');
end;

if length([event.value(ttl_idx)]) ~=LogDat.n
    error('number of events out of range');
end;

if any([event.samples(ttl_idx)] > nrSamples)
    error('events and CSC file out of range');
end;

params.event = event;

%%
switch params.mode
    case 'stimlocked'
        
        trl = [[event.samples(ttl_idx)]'-(2*params.pre*params.hdr.Fs) [event.samples(ttl_idx)]'+(2*params.pre*params.hdr.Fs)];
        trl(:,3) = -(2*params.pre*params.hdr.Fs)*ones(length(ttl_idx),1);
        
        del_idx = [find(sign(trl(:,1))==-1) find(trl(:,2)>event.samples(end)+2*params.post*params.hdr.Fs)];
        trl(del_idx,:) = [];
        LogDat.RT(del_idx) = [];
        %%
    case 'resplocked'
        
        trl = [[event.samples(ttl_idx2)]'-(2*params.pre*hdr.Fs) [event.samples(ttl_idx2)]'+(2*params.pre*hdr.Fs)];
        trl = double(trl);
        trl(:,3) = -(2*params.pre*hdr.Fs)*ones(length(ttl_idx2),1);
        
        del_idx = [find(sign(trl(:,1))==-1) find(trl(:,2)>event.samples(end)+2*params.post*hdr.Fs)];
        trl(del_idx,:) = [];
        LogDat.RT(del_idx) = [];
        ttl_idx = ttl_idx2;
end;

s = [event.samples(ttl_idx)]';

if max(s)>nrSamples
    error('timestamps and ttl events are out of range');
end;
%%
[params.TRL] = int64([s-params.pre*params.hdr.Fs s+params.post*params.hdr.Fs -params.pre*params.hdr.Fs*ones(length(s),1)]);
if max(max( params.TRL ))>nrSamples
    error('timestamps and ttl events are out of range');
end;
%% filter settings
Wp = [300 3000]./params.hdr.Fs;
[b,a] = butter(4,Wp);

params.Hd{1} = b;
params.Hd{2} = a;

%% open pool of workers
% if isempty(gcp('nocreate'))
%     %parpool(length(params.CSC_idx),'SpmdEnabled',false);
%     parpool(12,'SpmdEnabled',true);
% end;

%%
cd /home/rouxf/garbage/

%parfor ot = 1:length(params.CSC_idx)
for ot = 1:length(params.CSC_idx)
    %% select session to analyze
    
    params2 = params;
    params2.chanSel      = params.CSC_idx(ot);% range of channels to run analysis on
    [~,chan,~] = fileparts( params2.CSClabels{params2.chanSel} );
    
    params2.hdr.label = chan;
    params2.hdr.nChans = length(params2.hdr.label);
    
    chan
    
    [data,readme] = analyze_tuning_data2(params2);
    
    if ( params2.makeSpikes ==1 )
        chck = regexp(readme,'Spike');
        ix = [];
        for nt = 1:length(chck)
            ix(nt) = ~isempty(chck{nt});
        end;
        ix = find(ix);
        spike_data = data(ix);
        
        [outname] = [ params2.pID , '_' , 'TS' , num2str(params2.sesh) , '_spike_data_', chan , '_' , NLXts ,'_' , params2.mode, '.mat' ];
        
        par_save( [ params2.savepath , params.Nlxdat{params.sesh},  filesep, 'spike_dat' , filesep, outname ] , [] , {spike_data , params2 });
        
    end;
    
    if ( params2.makeLFP ==1 )
        chck = regexp(readme,'LFPdata');
        ix = [];
        for nt = 1:length(chck)
            ix(nt) = ~isempty(chck{nt});
        end;
        ix = find(ix);
        lfp_data = data(ix);
        
        [outname] = [ params2.pID , '_' , 'TS' , num2str(params2.sesh) , '_lfp_data_', chan , '_' , NLXts ,'_' , params2.mode, '.mat' ];
        
        par_save( [ params2.savepath , params.Nlxdat{params.sesh},  filesep, 'lfp_dat' , filesep, outname ] , [] , {lfp_data , params2 });
        
    end;
    
    if ( params2.makeMUA ==1 )
        chck = regexp(readme,'MUAdata');
        ix = [];
        for nt = 1:length(chck)
            ix(nt) = ~isempty(chck{nt});
        end;
        ix = find(ix);
        MUA_data = data(ix);
        
        [outname] = [ params2.pID , '_' , 'TS' , num2str(params2.sesh) , '_MUA_data_', chan , '_' , NLXts ,'_' , params2.mode, '.mat' ];
        
        par_save( [ params2.savepath , params.Nlxdat{params.sesh},  filesep, 'mua_dat' , filesep, outname ] , [] , {MUA_data , params2 });
        
    end;
    
end;
%%
delete('CSC_data_tmpCSC_*.dg_01.lab');
delete('CSC_data_tmpCSC_*.dg_01');
delete('spc_log.txt');

%%
%delete(gcp);
exit;