%%
cfg                     = [];
cfg.latency             = [-1 0];

[preDat] = ft_selectdata( cfg , EncDat );

cfg                     = [];
cfg.latency             = [0 1];

[postDat1] = ft_selectdata( cfg , EncDat );

cfg                     = [];
cfg.latency             = [2 3];

[postDat2] = ft_selectdata( cfg , EncDat );

dum1 = preDat;
dum2 = postDat1;
dum3 = postDat2;
for it = 1:length(preDat.trial)
    
    zP = zeros(size(preDat.trial{it}));
    zP = repmat(zP,[1 5]);
    
    dum1.trial{it} = [zP preDat.trial{it} zP];
    dum1.time{it} = 0:1/preDat.fsample:(length(dum1.trial{it})-1)/preDat.fsample;
    
    dum2.trial{it} = [zP postDat1.trial{it} zP];
    dum2.time{it} = 0:1/postDat1.fsample:(length(dum2.trial{it})-1)/postDat1.fsample;
    
    dum3.trial{it} = [zP postDat2.trial{it} zP];
    dum3.time{it} = 0:1/postDat2.fsample:(length(dum3.trial{it})-1)/postDat2.fsample;
    
end;

smth = 1/(length(dum1.time{1})/dum1.fsample);

cfg                     = [];
cfg.method              = 'mtmfft';
cfg.pad                 = 'maxperlen';
cfg.taper               = 'dpss';
cfg.tapsmofrq           = smth;

[powBaseL] = ft_freqanalysis( cfg , dum1 );
[powPostL1] = ft_freqanalysis( cfg , dum2 );
[powPostL2] = ft_freqanalysis( cfg , dum3 );

cfg                     = [];
cfg.frequency           = [1 15];

[powBaseL] = ft_selectdata( cfg, powBaseL );
[powPostL1] = ft_selectdata( cfg, powPostL1 );
[powPostL2] = ft_selectdata( cfg, powPostL2 );


dum1 = preDat;
dum2 = postDat1;
dum3 = postDat2;
for it = 1:length(preDat.trial)
    
    zP = zeros(size(preDat.trial{it}));
    zP = repmat(zP,[1 2]);
    
    dum1.trial{it} = [zP preDat.trial{it} zP];
    dum1.time{it} = 0:1/preDat.fsample:(length(dum1.trial{it})-1)/preDat.fsample;
    
    dum2.trial{it} = [zP postDat1.trial{it} zP];
    dum2.time{it} = 0:1/postDat1.fsample:(length(dum2.trial{it})-1)/postDat1.fsample;
    
    dum3.trial{it} = [zP postDat2.trial{it} zP];
    dum3.time{it} = 0:1/postDat2.fsample:(length(dum3.trial{it})-1)/postDat2.fsample;
    
end;

cfg                     = [];
cfg.method              = 'mtmfft';
cfg.pad                 = 'maxperlen';
cfg.taper               = 'dpss';
cfg.tapsmofrq           = 2.5;

[powBaseH] = ft_freqanalysis( cfg , dum1 );
[powPostH1] = ft_freqanalysis( cfg , dum2 );
[powPostH2] = ft_freqanalysis( cfg , dum3 );

cfg                     = [];
cfg.frequency           = [15 170];

[powBaseH] = ft_selectdata( cfg, powBaseH );
[powPostH1] = ft_selectdata( cfg, powPostH1 );
[powPostH2] = ft_selectdata( cfg, powPostH2 );

figure;
subplot(121);
plot(powPost1.freq,log(powPostL1.powspctrm)-log(powBaseL.powspctrm),'k');
axis tight;

subplot(122);
plot(powPostH1.freq,log(powPostH1.powspctrm)-log(powBaseH.powspctrm),'k');
axis tight;

cfg                     = [];
cfg.frequency           = [25 65];
cfg.avgoverfreq         = 'yes';

[dum1] = ft_selectdata( cfg , powBase2 );
[dum2] = ft_selectdata( cfg , powPost2 );

gPow = log(dum2.powspctrm)-log(dum1.powspctrm);
figure;
plot(sort(gPow),'ks-');

selIdx = find(gPow >0.8);

figure;
subplot(121);
plot(powPost1.freq,log(powPostL1.powspctrm(selIdx,:))-log(powBaseL.powspctrm(selIdx,:)),'k');
axis tight;
ylim([-.5 3]);

subplot(122);
plot(powPost2.freq,log(powPostH1.powspctrm(selIdx,:))-log(powBaseH.powspctrm(selIdx,:)),'k');
axis tight;
ylim([-.5 3]);