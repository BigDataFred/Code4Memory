epch1 = {};

%% get the event file data

p2d = '/media/rouxf/rds-share/Archive/MICRO/P05/TTL_NLX/2016-12-05_19-58-06/';

[dum] = getRawTTLs([p2d,'Events.nev']);
event.timestamp = dum(:,1);
event.value     = dum(:,2);

%% get onset events
picev = [event.value(:)];
sel_idx = max( find( picev==255 ) ):length( picev );
if isempty(sel_idx);
    sel_idx = 1:length( picev );
end;

ttl_idx = sel_idx(find(picev(sel_idx)==7));
ttl_idx(find([event.value(ttl_idx+1)]>7)) = [];% remove false bits
ttl_idx = ttl_idx';


%%
files = dir([p2d,'*.ncs']);
CSClabels = {};
for it = 1:length(files)
    
    [~,CSClabels{it},~] = fileparts(files(it).name);
    
end;

%% extract Fs and ADV
for ot = 1:length( CSClabels )
    
    FieldSelection(1) = 1;%timestamps
    FieldSelection(2) = 1;
    FieldSelection(3) = 1;%sample freq
    FieldSelection(4) = 1;
    FieldSelection(5) = 1;%samples
    ExtractHeader = 1;
    
    ExtractMode = 1; % 2 = extract record index range; 4 = extract timestamps range.
    %ModeArray(1)=[];
    %ModeArray(2)=[];
    
    [timestamps, ChanN, sampleFreq, nrADSamples, dataSamples , headerInfo] = Nlx2MatCSC_v3([ p2d, CSClabels{ot} ,'.ncs' ], FieldSelection, ExtractHeader, ExtractMode, []);
    dim = size(dataSamples);
    [idx] = getNlxHeaderParam(headerInfo,'ADBitVolts');% conversion factor
    scalef = str2double(headerInfo{idx}(12:end));
    
    %flatten
    dataSamples=dataSamples(:);
    dataSamples = dataSamples .* scalef .* 1e6;
    dataSamples = double(dataSamples);
    
    nrSamples = length(dataSamples);
    
    %%
    params.timestamps = timestamps;
    params.hdr.Fs = unique( sampleFreq );
    params.pre  = 0.0025;
    params.post = 0.0025;
    
    [trl] = convert_eventts2trl(event,params,ttl_idx);
    
    if max(max( trl ))>nrSamples
        error('timestamps and ttl events are out of range');
    end;
    
    %% sanity check
    
    for it = 1:size(trl,1)
        
        x = dataSamples(trl(it,1):trl(it,2)).*scalef.*1e6;
        epch1{ot}(it,1:length(x)) = x;
        
    end;
    
    if size(epch1{ot},2) ~= unique(diff(trl(:,1:2),[],2))+1;error('number of samples must fit epoch range');end;
    if size(epch1{ot},1) ~= size(trl,1);error('number of trial events is out of range');end;
    
end;

%%
epch2 = {};

%% get the event file data

p2d = '/media/rouxf/rds-share/Archive/MICRO/P05/TTL_NLX/2016-12-05_20-01-15/';

[dum] = getRawTTLs([p2d,'Events.nev']);
event.timestamp = dum(:,1);
event.value     = dum(:,2);

%% get onset events
picev = [event.value(:)];
sel_idx = max( find( picev==255 ) ):length( picev );
if isempty(sel_idx);
    sel_idx = 1:length( picev );
end;

ttl_idx = sel_idx(find(picev(sel_idx)==7));
ttl_idx(find([event.value(ttl_idx+1)]>7)) = [];% remove false bits
ttl_idx = ttl_idx';


%%
files = dir([p2d,'*.ncs']);
CSClabels = {};
for it = 1:length(files)
    
    [~,CSClabels{it},~] = fileparts(files(it).name);
    
end;

%% extract Fs and ADV
for ot = 1:length( CSClabels )
    
    FieldSelection(1) = 1;%timestamps
    FieldSelection(2) = 1;
    FieldSelection(3) = 1;%sample freq
    FieldSelection(4) = 1;
    FieldSelection(5) = 1;%samples
    ExtractHeader = 1;
    
    ExtractMode = 1; % 2 = extract record index range; 4 = extract timestamps range.
    %ModeArray(1)=[];
    %ModeArray(2)=[];
    
    [timestamps, ChanN, sampleFreq, nrADSamples, dataSamples , headerInfo] = Nlx2MatCSC_v3([ p2d, CSClabels{ot} ,'.ncs' ], FieldSelection, ExtractHeader, ExtractMode, []);
    dim = size(dataSamples);
    [idx] = getNlxHeaderParam(headerInfo,'ADBitVolts');% conversion factor
    scalef = str2double(headerInfo{idx}(12:end));
    
    %flatten
    dataSamples=dataSamples(:);
    dataSamples = dataSamples .* scalef .* 1e6;
    dataSamples = double(dataSamples);
    
    nrSamples = length(dataSamples);
    
    %%
    params.timestamps = timestamps;
    params.hdr.Fs = unique( sampleFreq );
    params.pre  = 0.0025;
    params.post = 0.0025;
    
    [trl] = convert_eventts2trl(event,params,ttl_idx);
    
    if max(max( trl ))>nrSamples
        error('timestamps and ttl events are out of range');
    end;
    
    %% sanity check
    
    for it = 1:size(trl,1)
        
        x = dataSamples(trl(it,1):trl(it,2)).*scalef.*1e6;
        epch2{ot}(it,1:length(x)) = x;
        
    end;
    
    if size(epch2{ot},2) ~= unique(diff(trl(:,1:2),[],2))+1;error('number of samples must fit epoch range');end;
    if size(epch2{ot},1) ~= size(trl,1);error('number of trial events is out of range');end;
    
end;

%% sanity check
figure;
for it = 1:length(epch1)      
    subplot(6,8,it);
    hold on;
    imagesc(linspace(-params.pre,params.post,size(epch1{it},2)),1:size(epch1{it},1),abs(epch1{it}));
    %plot(linspace(-params.pre,params.post,size(epch1{it},2)),mean(epch1{it},1),'b');
    %plot(linspace(-params.pre,params.post,size(epch1{it},2)),mean(epch1{it},1),'r.');
    axis tight;
end;
%% sanity check
%figure;
%subplot(4,1,1:3);
%imagesc(linspace(-params.pre,params.post,size(epch2{ot},2)),1:size(epch2{ot},1),epch2{ot});
%subplot(4,1,4);
%plot(linspace(-params.pre,params.post,size(epch2{ot},2)),mean(epch2{ot},1));