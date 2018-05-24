function [delIx] = detectIEDevents( lfpDat )

%% band-pass filter
cfg                     = [];
cfg.bpfilter            = 'yes';
cfg.bpfreq              = [25 80];

[bfDat] = ft_preprocessing( cfg , lfpDat );

%% extract filtered envelope
cfg                     = [];
cfg.hilbert             = 'abs';

[bfDat] = ft_preprocessing( cfg , bfDat );

%%
delIx = [];
for it = 1:length(bfDat.trial)
    
    %amplitude envelope filtered
    x = bfDat.trial{it}.^2;
    m = median(x,2);
    m = m*ones(1,size(x,2));        
    sd = iqr(x')';
    sd = sd*ones(1,size(x,2));
    [z1] = x>m+2.*sd;%(x-m)./sd;
    
    % slope envelope filtered data
    x = gradient(bfDat.trial{it}).^2;
    m = median(x,2);
    m = m*ones(1,size(x,2));
    sd = iqr(x')';
    sd = sd*ones(1,size(x,2));
    [z2] = x>m+2.*sd;%(x-m)./sd;
    
    % slope raw data
    x = gradient(lfpDat.trial{it}).^2;
    m = median(x,2);
    m = m*ones(1,size(x,2));        
    sd = iqr(x')';
    sd = sd*ones(1,size(x,2));
    [z3] = x>m+2.*sd;%(x-m)./sd;
     
    % amplitude raw data
    x = abs(hilbert(lfpDat.trial{it})).^2;
    m = median(x,2);
    m = m*ones(1,size(x,2));        
    sd = iqr(x')';
    sd = sd*ones(1,size(x,2));
    [z4] = x>m+2.*sd;%(x-m)./sd;       
    
    thrsh = z1+z2+z3+z4;
    thrsh = thrsh > 3;
    
%     ix = find(thrsh);
%     [chan,smp] = ind2sub(size(thrsh),ix);
%     [~,sIdx] = sort(chan);
%     chan = chan(sIdx);
%     smp = smp(sIdx);
%     
%     chanID = unique(chan);
%     sel = [];
%     for jt = 1:length(chanID)        
%         ix2 = find(chan == chanID(jt));
%         dt = diff([smp(ix2)])';
%         f = find(diff([false, dt==1, false])~=0);
%         out = f(2:2:length(f))-f(1:2:length(f));
%         thr = (out/lfpDat.fsample > 0.02);% & (out/lfpDat.fsample < 0.1);
%         if ( any(thr~=0) )           
%             sel = [sel jt];
%         end;
%     end;
    
    if max(max(thrsh)) ~=0%~isempty(sel) %&& (max(max(orig)) - min(min(orig)) >8)
        delIx = [delIx it];
    end;
    
end;
