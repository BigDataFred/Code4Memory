function [trl] = convert_eventts2trl(event,params,ttl_idx)

%% convert the event timestamps to ms
x1 = [event.timestamp] ;
x1 = x1./1e6;% convert from micro-s to s

if ~isempty(ttl_idx)
    x1 = x1(ttl_idx);% keep only those events that match requirements
end;

%FIXME
% if strcmp(params.mode,'resplocked')
%     x1 = x1 + params.RTs+3;% RTs must be in s
% end;

%% interpolate the data timestamps
x2 = params.timestamps;
[tsi] = timeStampinter(x2);%interpolate timestamps
x2 = tsi./1e6;% convert from micro-s to s

%% get the epoch time range
dt = zeros(length(x1),2);
for it = 1:length(x1)
    dt(it,:) = [x1(it)-params.pre x1(it)+params.post];%duration of trial including baseline
end;

%% find the dataSamples that correspond to time range of each epoch
ix = cell(1,length(x1));
n  = zeros(1,length(x1));
nsamp = unique( single( params.hdr.Fs * unique(diff(dt,[],2)) ) );
for it = 1:size(dt,1)
    
    ix{it} = find( x2 >= dt(it,1) & x2 <= dt(it,2) );
    if length(ix{it}) > nsamp
        ix{it} = ix{it}(1:nsamp);%make sure to keep equal number of samples
    end;
    n(it)  = length( ix{it} );
    
end;

if length(unique(n)) ~= 1%quick sanity check
    error('number of samples is out of range');
end;

%% build the trl matrix
trl = zeros(length(ix),3);
for it = 1:length(ix)
    
    trl(it,:) = [min(ix{it}) max(ix{it}) -params.pre*params.hdr.Fs];
    
end;

