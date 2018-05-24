function [trl_dat] = make_spike_trials(params, spike)

u = spike.unit{1};
[d1,d2] = size(u);
if d1>d2
    u = u';
end;

ts = double(spike.timestamp{1}) -  double(params.hdr.FirstTimeStamp);
ts2 = double(spike.timestamp{1});% -  double(params.hdr.FirstTimeStamp);

ttl_idx = [params.event(:).value];% all TTLs
ttl_idx = find(ttl_idx ==7);% get indexes of events of interest
picon = [params.event(:).timestamp];
picon = picon -  double( params.hdr.FirstTimeStamp );
picon = picon(ttl_idx);
waveform = [];

%%
trial     = []; time      = []; timestamp = []; trialtime = []; unit = [];  waveform = [];

%%
for it = 1:length(picon)
    
    x = ts - picon(it);
    x = double(x).*1e-6;
    
    sel = [];
    sel = find(x>=-params.pre & x <=params.post);
    
    if ~isempty(sel)
        trial     = [trial it*ones(1,length(sel))];
        time      = [time x(sel)];
        timestamp = [timestamp ts2(sel)];
        unit      = [unit u(sel)];
        
        if length(sel)>1
            waveform = [waveform;squeeze(spike.waveform{1}(:,:,sel))'];
        else
            waveform = [waveform;squeeze(spike.waveform{1}(:,:,sel))];
        end;
    end;
    sel = [];
end;
waveform = waveform';

trialtime = [-params.pre*ones(length(picon),1) params.post*ones(length(picon),1)];

%%
trl_dat.trial       = trial;
trl_dat.time        = time;
trl_dat.timestamp   = timestamp;
trl_dat.trialtime   = trialtime;
trl_dat.unit        = unit;
trl_dat.waveform    = waveform;
