function [phi,pow] = extractPhaseLFP(LFPdat,foi,Fs,flag)
%%
phi = {};
pow = {};

%%
if isstruct(LFPdat)
    %case data is in ft format
else
    % switch format from Chronux to fieldtrip
    dum = [];
    dum.trial = cell(1,size(LFPdat{1},2));
    dum.time = cell(1,size(LFPdat{1},2));
    dum.label = cell(1,length(LFPdat));
    dum.Fs = Fs;
    Nsamp = size(LFPdat{1},1);
    for it = 1:length(LFPdat)% loop over the channels
        for jt = 1:size(LFPdat{it},2)%loop over the trials
            dum.trial{jt}(it,:) = [fliplr(LFPdat{it}(:,jt)') LFPdat{it}(:,jt)' fliplr(LFPdat{it}(:,jt)')];
            dum.time{jt} = 0:1/Fs:(Nsamp*3-1)/Fs;
        end;
        dum.label(it) = {['dummyChan',num2str(it)]};
    end;
    LFPdat = dum;
end;

if flag(1) == 1
    cfg                     = [];
    cfg.method              = 'wavelet';
    cfg.foi                 = foi;
    cfg.toi                 = dum.time{1};
    cfg.width               = 4;
    cfg.output              = 'fourier';
    cfg.pad                 = 'nextpow2';
    
    [phi] = ft_freqanalysis( cfg , LFPdat );
    
    phi.fourierspctrm = phi.fourierspctrm(:,:,:,Nsamp*2+1:Nsamp*3);
    %chck = squeeze(phi.fourierspctrm(1,1,1,:));
   
    %if any(isnan(angle(chck)))
    %    error('phase estimate contains NaNs');
    %end;
    x = cell(1,length(phi.label));
    for it = 1:length(phi.label)
        x{it} = angle(squeeze(phi.fourierspctrm(:,it,:,:)));
    end;
    phi = x;
end;

if flag(2) == 1
    cfg                     = [];
    cfg.method              = 'wavelet';
    cfg.foi                 = foi;
    cfg.toi                 = dum.time{1};
    cfg.width               = 4;
    cfg.output              = 'pow';
    cfg.pad                 = 'nextpow2';
    
    [pow] = ft_freqanalysis( cfg , LFPdat );
    
    %pow.powspctrm = pow.powspctrm(:,:,Nsamp*3+1:Nsamp*4);
    %chck = squeeze(pow.powspctrm(1,1,1));
    %if any(isnan(chck))
    %    error('power estimate contains NaNs');
    %end;        
    
    x = cell(1,length(pow.label));
    for it = 1:length(pow.label)
        x{it} = squeeze(pow.powspctrm(it,:,:));
    end;
    pow = x;
end;

return;
