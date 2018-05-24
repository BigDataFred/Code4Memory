%% add fieldtrip to matlab path
addpath('/media/rouxf/rds-share/Common/fieldtrip-20170115/');
ft_defaults;

%% reat the triggers from the macro channels
path2MAC = '/media/rouxf/rds-share/Archive/MICRO/P07/visualStim/';
macFN = 'IG_gratings.edf';

cfg                     = [];
cfg.dataset             = [path2MAC,macFN];
cfg.channel             = {'TRIG'};

[trigDat] = ft_preprocessing(cfg);% read the trigger channel from file

%% find samples corresponding to triggers
dt = diff(trigDat.trial{1});% difference between adjacent samples
[ixON] = find(sign(dt)==-1)+1;% onsets of trigger
[ixOFF] = find(sign(dt)==1);% offsets of trigger
ixOFF = ixOFF(1:length(ixON));
  
onSamp =ixON';
preSamp = 1*trigDat.fsample;% samples to substract to get the baseline
postSamp = 3*trigDat.fsample;% samples to add to get the end of the poststim window


%% segment the macro data into epochs around the trigger events
cfg                     =[];
cfg.dataset             = [path2MAC,macFN];
cfg.channel             = {'all', '-D*', '-Event', '-TRIG', '-OSAT', '-PR', '-Pleth'};
cfg.lpfilter            = 'yes';
cfg.lpfilttype          = 'fir';
cfg.lpfreq              = 200;
cfg.detrend             = 'yes';
cfg.demean              = 'yes';
cfg.continuous          = 'yes';
cfg.trl                 = round([onSamp-preSamp onSamp+postSamp -preSamp*ones(length(onSamp),1)]);

[macroDat] = ft_preprocessing( cfg );

%% extract the channels of interest for subsequent analysis
cfg                     = [];
cfg.trials              = 1:95;
cfg.channel             = {'R C_*'};

[chanSel] = ft_selectdata( cfg , macroDat );

%% browse the data manually and mark trials with artefacts
cfg                     = [];
cfg.viewmode            = 'vertical';

[cfg] = ft_databrowser( cfg , chanSel );

%% reject the trials that were marked by hand
[chanSel] = ft_rejectartifact( cfg , chanSel );

%% common average re-referecing
x = cell(1,length(chanSel.trial));

for it = 1:length(chanSel.trial)% loop over trials
    
    dat = chanSel.trial{it};
    
    for jt = 1:length(chanSel.label)%loop over channels
        d = dat(jt,:) - mean(dat,1);% substract the average across all channels
        x{it}(jt,:) = d;
    end;
    
end;

% create dummy structure for subsequent analysis
dum = [];
dum.time = chanSel.time;
dum.trial = x;
dum.label = chanSel.label;

%% apply low-pass filter and compute ERP
cfg                     = [];
cfg.lpfilter            = 'yes';
cfg.lpfreq              = [30];

[dum2] = ft_preprocessing( cfg , dum );


itc = [];
for it = 1:length(dum2.trial)
    for jt = 1:length(dum2.label)
    
        itc(it,jt,:) = hilbert( dum2.trial{it}(jt,:));
        
    end;
end;

itc = itc./abs(itc);
itc = squeeze(1/size(itc,1)*abs(sum(itc,1)));

cfg                     = [];
cfg.keeptrials          = 'yes';

[tlck] = ft_timelockanalysis( cfg , dum2 );

%% make single trial ERP-image
n =length(tlck.label)/8;

sel = 1:size(tlck.trial,1);
figure;
for it = 1:length(tlck.label)
    subplot(2,4,it);
    a(it) =gca;
    imagesc(tlck.time,1:size(tlck.trial,1),squeeze(tlck.trial(sel,it,:)));   
    caxis([min(min(squeeze(tlck.trial(sel,it,:))))/3 max(max(squeeze(tlck.trial(sel,it,:))))/3]);
    axis xy;
    xlim([-.5 2.5]);
    x = tlck.label{it};
    x(regexp(x,'_')) = [];
    x(regexp(x,' ')) = [];
    title(x);
    hold on;
    plot([0 0],[1 size(tlck.trial,1)],'w');
    plot([1 1],[1 size(tlck.trial,1)],'w');
    plot([2 2],[1 size(tlck.trial,1)],'w');
    x = itc(it,:);
    %x = (x-min(x))./(max(x)-min(x));
    x = 50+(25.*x);
    plot(tlck.time,x,'Color',[.9 0 0]);
    xlim([-.5 2.5]);
    
end;

for it = 1:length(a)
    xlabel(a(it),'Time (s)');
    ylabel(a(it),'Trial #');
end;

%% do bipolar referencing
x = cell(1,length(chanSel.trial));

idx =[];
c = 0;
for it = 1:2:length(chanSel.label)
    c = c+1;
    idx(c,:) = [it it+1];% indexes of adjacent channels
end;


for it = 1:length(chanSel.trial)
    
    dat = chanSel.trial{it};
    
    for jt = 1:size(idx,1)        
        d = dat(idx(jt,1),:) - dat(idx(jt,2),:);% bipolar reference
        x{it}(jt,:) = d;
    end;
    
end;

l = cell(1,size(idx,1));
for jt = 1:size(idx,1)
    l(jt) = {['bipolar_',num2str(jt)]};
end;

% create dummy structure for subsequent analysis
dum = [];
dum.time = chanSel.time;
dum.trial = x;
dum.label = l;

% %%
% for it = 1:length(dum.trial)
%     for jt = 1:length(dum.label)
%         dum.trial{it}(jt,:) = dum.trial{it}(jt,:) - tlck.avg(jt,:);
%     end;
% end;

%% do the time-frequency analysis
cfg                         = [];
cfg.method                  = 'mtmconvol';
cfg.pad                     = 'maxperlen';
cfg.taper                   = 'dpss';
cfg.toi                     = chanSel.time{1}:0.05:chanSel.time{1}(end);
cfg.foi                     = [1:19 20:4:200];
cfg.tapsmofrq(cfg.foi<20)   = ones(1,length(find(cfg.foi<20)));
cfg.tapsmofrq(cfg.foi>=20)  = 10*ones(1,length(find(cfg.foi>=20)));
cfg.t_ftimwin(cfg.foi<20)   = ones(1,length(find(cfg.foi<20)));
cfg.t_ftimwin(cfg.foi>=20)  = 0.25*ones(1,length(find(cfg.foi>=20)));
cfg.keeptrials              = 'yes';

[pow] = ft_freqanalysis( cfg , dum );

%% compute baseline correction
cfg                     = [];
cfg.baseline            = [-.5 -.1];
cfg.baselinetype        = 'relchange';

[bcpow]  = ft_freqbaseline( cfg , pow );

%% visualise time-frequency maps
n =length(pow.label)/4;

sel = 1:size(pow.powspctrm,1);
figure;
a = [];
for it = 1:length(pow.label)
    subplot(n,4,it);
    a(it) =gca;
    Y = squeeze(mean(mean((bcpow.powspctrm(sel,it,:,:)),2),1));
    pcolor(bcpow.time,bcpow.freq,Y);shading interp;
    %caxis([min(min(Y))/10 max(max(Y))/10]);
    hold on;
    plot([0 0],[bcpow.freq(1) bcpow.freq(end)],'w');
    plot([1 1],[bcpow.freq(1) bcpow.freq(end)],'w');
    plot([2 2],[bcpow.freq(1) bcpow.freq(end)],'w');
    axis xy;axis tight;
    xlim([-.5 2.5]);ylim([20 200]);
    caxis([-1 1]);
    x = pow.label{it};
    x(regexp(x,'_')) = [];
    x(regexp(x,' ')) = [];
    title(x);
end;

for it = 1:length(a)
    xlabel(a(it),'Time (s)');
    ylabel(a(it),'Frequency (Hz)');
end;
