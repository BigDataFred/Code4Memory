function [PPC,freqAx,pval] = computePPC(data_lfp,data_spk,Fs,timeGrid,bootStrap)
%%
ntrl = size(data_lfp,2);

[dum] = plugDat2FT(data_lfp,data_spk,timeGrid,Fs);

cfg1                     = [];
cfg1.method              = 'mtmfft'; 
cfg1.timwin              = [-.5 .5];
cfg1.spikechannel        = 'SpkChan1';
cfg1.channel             = 'LFPchan1';
cfg1.foilim              = [1 30];
cfg1.taper               = 'dpss';
cfg1.tapsmofrq           = 5;

[sts] = ft_spiketriggeredspectrum( cfg1, dum);

cfg2                         = [];
cfg2.method                  = 'ppc0';
cfg2.spikechannel             = 'SpkChan1';
cfg2.channel                 = 'LFPchan1';
cfg2.avgoverchan             = 'unweighted';
cfg2.timwin                  = 'all';

[statSts] = ft_spiketriggeredspectrum_stat( cfg2, sts);
PPC = eval(['statSts.',cfg2.method]);
freqAx = statSts.freq;

niter = 500;
nfreq = length(statSts.freq);
randPPC = zeros(niter,nfreq);
pval = zeros(1,nfreq);
if strcmp(bootStrap,'y')
    parfor it = 1:niter
        
        randTrl = randperm(ntrl);
        randDat  = data_lfp(:,randTrl);
        
        [dum] = plugDat2FT(randDat,data_spk,timeGrid,Fs);
        
        [sts] = ft_spiketriggeredspectrum( cfg1, dum);
        
        [statSts] = ft_spiketriggeredspectrum_stat( cfg2, sts);
        randPPC(it,:) = statSts.ppc0;
        
    end;
        
    parfor it = 1:nfreq
        pval(it) = length(find(randPPC(:,it) >= PPC(it)))/niter
    end;
end;
return;

function [dum] = plugDat2FT(data_lfp,data_spk,timeGrid,Fs)
x = data_lfp;
y = data_spk;
ntrl = size(data_lfp,2);

dum = [];
dum.fsample = Fs;
dum.time = {};
dum.trial = {};
dum.label = {'LFPchan1','SpkChan1'};

for it = 1:ntrl
    dum.time{it} = timeGrid;
    dum.trial{it} = [x(:,it)';y(:,it)'];
end;

