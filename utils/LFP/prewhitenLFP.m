function [pwDat] = prewhitenLFP(LFPdat,method,Fs)

if isstruct(LFPdat)
    f =1;
else
    f=2;
    if strcmp(method,'ARMfilt')
        dum                     = [];
        dum.fsample             = Fs;
        dum.time                = {};
        dum.trial               = {};
        dum.label               = {};
        for it = 1:length(LFPdat)
            dum.label(it) = {['dummyChan',num2str(it)]};
            for jt = 1:size(LFPdat{it},2)
                dum.trial{jt}(it,:) = LFPdat{it}(:,jt);
                nsamp = size( LFPdat{it},1 );
                dum.time{jt}        = 0:1/Fs:(nsamp-1)/Fs;
            end;
        end;
        LFPdat = dum;
    end;
end;

switch method
    
    case 'dx/dt' % take the derivative
        
        if f ==1
            for it = 1:length(LFPdat.trial)
                [Fx,~] = gradient(LFPdat.trial{it}); % Fx = horizontal derivative (time), Fy = vertical derivative (channels)
                LFPdat.trial{it} = Fx;% only keep derivative over time
            end;
        elseif f==2
            for it = 1:length(LFPdat)
                [~,Fy] = gradient(LFPdat{it});% Fx = horizontal derivative (channels), Fy = vertical derivative (time)
                LFPdat{it} = Fy;
            end;
        end;
        
    case 'ARMfilt' % estimate coefficients of ARM model and use with filter
        %selIdx = find(LFPdat.time{1} >=-1 & LFPdat.time{1} <0);
        selIdx = 1:length(LFPdat.time{1});
        
        ntrl =  length(LFPdat.time);
        nsamp = length(selIdx);
        idx = 1:nsamp;
        concat = zeros(length(LFPdat.label),nsamp*ntrl);
        for it = 1:length(LFPdat.trial)
            concat(:,idx) = LFPdat.trial{it}(:,selIdx);
            idx = idx + nsamp;
        end;
        
        dt = LFPdat.time{1}(2)-LFPdat.time{1}(1);
        Fs = 1/dt;
        
        dum                     = LFPdat;
        dum.trial               = {};
        dum.time                = {};
        dum.trial{1}            = concat;
        dum.time{1}             = 0:1/Fs:(length(concat)-1)/Fs;
        
        cfg                     = [];
        cfg.order               = 2;
        cfg.toolbox             = 'bsmart';
        cfg2                    = [];
        cfg2.method             = 'mvar';
        
        A = [];
        for jt = 1:length(dum.label)
            
            cfg.channel         = dum.label(jt);
            [mdata] = ft_mvaranalysis(cfg,dum);
            
            [mfreq] = ft_freqanalysis(cfg2,mdata);
            
            A(jt,:) = [1 -squeeze(mdata.coeffs(1,1,:))'];
        end;
        
        %% filter LFP with autoregressive model
        for jt = 1:length(LFPdat.label)
            for it = 1:length(LFPdat.trial)
                
                y = filtfilt(A(jt,:),1,LFPdat.trial{it}(jt,:));
                if any(isnan(y))
                    error('toto');
                else
                    LFPdat.trial{it}(jt,:) = y;
                end;
                
            end;
        end;
end;
if f==1 || ( (f==2) && (strcmp(method,'dx/dt')) )
    pwDat = LFPdat;        
elseif f==2 && strcmp(method,'ARMfilt')
    pwDat = cell(1,length(LFPdat.label));
    for it = 1:length(LFPdat.label)
        pwDat{it} = zeros(size(LFPdat.trial{1}'));
        for jt = 1:length(LFPdat.trial)
            pwDat{it}(:,jt) = LFPdat.trial{jt}(it,:);
        end;
    end;
end;