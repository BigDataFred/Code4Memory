function [data_lfp,data_lfp_interp,data_spk,data_spk2] = prepSPKandLFPdat4Chronux(spkDat,MWdat,lfpMWselIx,spkMWselIx,trlIdx,toi,Fs)

data_spk        = [];
data_spk2       = [];
data_lfp        = [];
data_lfp_interp = [];

tlck = [];
if ~isempty(lfpMWselIx)
    cfg                     = [];
    cfg.keeptrials          = 'yes';
    cfg.channel             = MWdat.label(lfpMWselIx);
    cfg.preproc.demean       = 'yes';
    cfg.preproc.detrend      = 'yes';
    
    [tlck] = ft_timelockanalysis( cfg ,MWdat );
    
    cfg                     =[];
    cfg.latency             = [toi(1) toi(2)];
    
    [tlck] = ft_selectdata( cfg, tlck );
    data_lfp = squeeze(tlck.trial)';
    dt = [min(tlck.time) max(tlck.time)];
    dt = dt.*1e3;
    dt = dt(1):dt(2);
    [data_lfp_interp] = zeros(length(dt),length(trlIdx));
end;

if ~isempty(spkMWselIx)
                
    trl = spkDat{spkMWselIx}.trial{1};
    ts = spkDat{spkMWselIx}.time{1};
    
    if length(trl) ~= length(ts)
        error('trial and timestamp assignment must be equal');
    end;
    
    trl(ts<toi(1) & ts>toi(end)) = [];
    ts(ts<toi(1) & ts>toi(end)) = [];
    ts = ts.*1e3;
    if length(trl) ~= length(ts)
        error('trial and timestamp assignment must be equal');
    end;
    
    dt = toi.*1e3;
    dt = dt(1):dt(2);
    [data_spk]        = zeros(length(dt),length(trlIdx));
    [data_spk2]       = struct;
    for it = 1:length(trlIdx)
        x = ts(trl == trlIdx(it));
        x(x<min(dt)) =[];
        x(x > max(dt)) = [];
        data_spk(:,it) = hist(x,dt);
        data_spk2(it).times = x./1e3;
        if ~isempty(lfpMWselIx)
            [data_lfp_interp(:,it)] = interpLFP(data_lfp(:,it),data_spk(:,it),[2 6].*1e-3,Fs,'linear');
        end;
    end;    

end;