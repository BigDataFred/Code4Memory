%% reset workspace
clear;
clc;

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
params.CSC_idx = [ 1:24 ];
params.pID = 'P04';
params.sesh = [ 4 ];
params.p2Nlxdata = ['/media/rouxf/rds-share/Archive/MICRO/',params.pID,'/fvSpEM/'];%path 2 Nlx data
params.savepath = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',params.pID,'/fvSpEM/'];%path for saving
params.rpath = '/home/rouxf/';

params.mode         = 'stimlocked';%'resplocked';%
params.pre          = 3;% baseline time range
params.post         = 3;% params.post stim time range

% set flags
params.makeLFP      = 1;
params.makeSpikes   = 1;
params.makeMUA      = 1;

chck = dir( params.savepath );
if isempty( chck )
    mkdir( params.savepath )
end;

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

%%
ix      = regexp(params.p2Nlxdata,'/');
NLXts     = params.p2Nlxdata(ix(end-1):end);
NLXts(regexp(NLXts,'/')) = [];

%% get the logfile data
[sdat] = create_session_dat(params.pID,params.sesh);

params1.p = sdat.p2lf;%path
params1.fn = sdat.lfn;%filename
params1.ntrl = 49;%number of trials
params1.ncol = 9;% number of columns in logfile

%same as above
params2.p = sdat.p2lf;
params2.fn = sdat.lfn;
params2.ntrl = 49;
params2.ncol = 12;

%read logfile for Encoding
[LogDat1] = getNewLogDataEM( params1, 'ENC' );
%read logfile for Retrieval
[LogDat2] = getNewLogDataEM( params2, 'RET' );

[RTs] = str2double( LogDat1.log(:,end) );%get RTs
params.RTs = RTs;

%find hits, misses and stim categories
cat = LogDat1.log(:,3:4);
c = {};
for jt = 1:size( cat,1 )
    
    c{jt} = [ cat{jt,1}(1) cat{jt,2}(1) ];% reads out info from logdata
    
end;

ix    = {};
ix{1} = find( strcmp(c,'fp') );% indices corresponding to face-place tr
ix{2} = find( strcmp(c,'pp') );% indices corresponding to place-place tr
ix{3} = find( strcmp(c,'ff') );% indices corresponding to face-face tr

ix{4} = find( sum( [str2double(LogDat2.log(:,5:6))] ,2) ==2);% both images were correctly remebered
ix{5} = find( sum( [str2double(LogDat2.log(:,5:6))] ,2) ==1);% only 1/2 images were correctly remebered
ix{6} = find( sum( [str2double(LogDat2.log(:,5:6))] ,2) ==0);%no images were correctly remembered

ix_readme = { 'face-place' 'place-place' 'face-face' 'correct-both' 'miss-one' 'miss-both' };% save info for later
LogDat_readme = { 'LogDat1:Enc', 'LogDat2:Ret' };

%%
[~,fn,~] = fileparts(sdat.lfn);
outname = [ fn , '_LogDat.mat' ];

save( [ params.savepath , params.Nlxdat{params.sesh},  filesep, 'log_dat' , filesep, outname ] , 'LogDat1','LogDat2','ix','ix_readme','LogDat_readme','RTs' );

%% get the event file data
[dum] = getRawTTLs([params.p2Nlxdata,'Events.nev']);
event.timestamp = dum(:,1);
event.value     = dum(:,2);
%[event] = ft_read_event([params.p2Nlxdata,'Events.nev']);

%% get onset events
picev = [event.value(:)];
sel_idx = max( find( picev==255 ) ):length( picev );

if isempty(sel_idx)
    sel_idx = 1:length(picev);
end;

ttl_idx = sel_idx(find(picev(sel_idx)==7));
ttl_idx(find([event.value(ttl_idx+1)]>7)) = [];% remove false bits
ttl_idx = ttl_idx';

if sum( [ length(LogDat1.log) length(LogDat2.log) ] ) ~= length(ttl_idx)
    error('number of events is out of range');
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

[idx] = getNlxHeaderParam(headerInfo,'AcqEntName');% CSC name
params.hdr.label = headerInfo{idx}(13:end);
params.hdr.nChans = length(params.hdr.label);

params.hdr.FirstTimeStamp = int64(params.timestamps(1));% first ts
params.hdr.nSamples = nrSamples;
params.hdr.TimeStampPerSample = (params.timestamps(2)-params.timestamps(1))/512;

%% convert event ts to samples
x = [event.timestamp(:)];
n = length( x );
for it = 1:n
    event.samples(it) = (x(it) - double(params.hdr.FirstTimeStamp))./params.hdr.TimeStampPerSample + 1;    
end;
event.samples = event.samples';

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

if length([event.value(ttl_idx)]) ~= sum( [ length(LogDat1.log) length(LogDat2.log) ] )
    error('number of events out of range');
end;

if any([event.samples(ttl_idx)] > nrSamples)
    error('events and CSC file out of range');
end;

params.event = event;

%%
% switch params.mode
%     case 'stimlocked'
%         
%         trl = [[event.samples(ttl_idx)]-(2*params.pre*params.hdr.Fs) [event.samples(ttl_idx)]+(2*params.post*params.hdr.Fs)];
%         trl(:,3) = -(2*params.pre*params.hdr.Fs)*ones(length(ttl_idx),1);
%         
%         del_idx = [find(sign(trl(:,1))==-1) find(trl(:,2)>event.samples(end)+2*params.post*params.hdr.Fs)];
%         trl(del_idx,:) = [];
%         RTs(del_idx) = [];
%         %%
%     case 'resplocked'
%         
%         trl = [[event.samples(ttl_idx2)]-(2*params.pre*hdr.Fs) [event.samples(ttl_idx2)]+(2*params.post*hdr.Fs)];
%         trl = double(trl);
%         trl(:,3) = -(2*params.pre*hdr.Fs)*ones(length(ttl_idx2),1);
%         
%         del_idx = [find(sign(trl(:,1))==-1) find(trl(:,2)>event.samples(end)+2*params.post*hdr.Fs)];
%         trl(del_idx,:) = [];
%         RTs(del_idx) = [];
%         ttl_idx = ttl_idx2;
% end;
% 
% s = [event.samples(ttl_idx)];
% 
% if max(s)>nrSamples
%     error('timestamps and ttl events are out of range');
% end;

%%
[trl] = convert_eventts2trl(event,params,ttl_idx);

params.TRL = trl;
if max(max( params.TRL ))>nrSamples
    error('timestamps and ttl events are out of range');
end;
%% sanity check
epch = [];

for it = 1:size(trl,1)
    
    x = dataSamples(trl(it,1):trl(it,2)).*params.hdr.scalef.*1e6;
    epch(it,1:length(x)) = x;
    
end;

if size(epch,2) ~= unique(diff(trl(:,1:2),[],2))+1;error('number of samples must fit epoch range');end;
if size(epch,1) ~= size(trl,1);error('number of trial events is out of range');end;

%% sanity check
figure;
subplot(4,1,1:3);
imagesc(linspace(-params.pre,params.post,size(epch,2)),1:size(epch,1),epch);
subplot(4,1,4);
plot(linspace(-params.pre,params.post,size(epch,2)),mean(epch,1));

%%
% [params.TRL] = int64([s-params.pre*params.hdr.Fs s+params.post*params.hdr.Fs -params.pre*params.hdr.Fs*ones(length(s),1)]);
% if max(max( params.TRL ))>nrSamples
%     error('timestamps and ttl events are out of range');
% end;

%% filter settings
Wp = [300 3000]./params.hdr.Fs;
[b,a] = butter(4,Wp);

params.Hd{1} = b;
params.Hd{2} = a;

%% open pool of workers
if isempty(gcp('nocreate'))
    parpool(length(params.CSC_idx)/3,'SpmdEnabled',false);
end;

%%
cd( [ params.rpath,'garbage' ]);

parfor ot = 1:length(params.CSC_idx)
%for ot = 1:length(params.CSC_idx)
    %% select session to analyze
    
    params2 = params;
    params2.chanSel      = params.CSC_idx(ot);% range of channels to run analysis on
    [~,chan,~] = fileparts( params2.CSClabels{params2.chanSel} );
    chan
    
    [data,readme] = analyze_EM_data(params2);
    
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
    
    if ( params2.makeSpikes ==1 )
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
    
    if ( params2.makeSpikes ==1 )
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
delete(gcp);
exit;



