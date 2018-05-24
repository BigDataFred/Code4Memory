%%
clear all;
clc;
%%
if isempty(gcp('nocreate'))
    parpool(32);
end;
%% set path defs
restoredefaultpath;
addpath('~rouxf/fieldtrip-20161009/');
addpath(genpath('~rouxf/AnalysisFred/'));
ft_defaults;
%%
sesh = 1;
pID = 'P02';

TS{1} = '2016-07-09_11-13-28';
LF{1} = 'P02_TS01_log_ctune_09-Jul-2016_11_40_xx.txt';
 
% TS{1} = '2016-10-16_15-15-03';
% TS{2} = '2016-10-17_14-04-16';
% TS{3} = '2016-10-17_16-37-42';
% TS{4} = '2016-10-18_16-07-23';
% TS{5} = '2016-10-19_14-34-32';
% TS{6} = '2016-10-20_14-26-32';
% 
% LF{1} = 'P04_TS01_TS01_log_ctune_16102016_15_43_27.txt';
% LF{2} = 'P04_TS02_TS02_log_ctune_17102016_14_4_32.txt';
% LF{3} = 'P04_TS03_TS03_log_ctune_17102016_17_15_25.txt';
% LF{4} = 'P04_TS04_TS04_log_ctune_18102016_16_6_3.txt';
% LF{5} = 'P04_TS05_TS05_log_ctune_19102016_15_35_24.txt';
% LF{6} = 'P04_TS06_TS06_log_ctune_20102016_14_25_16.txt';
%% files with timestamps
p2f = ['/home/rouxf/res/tuning/',pID,'/',TS{sesh},'/'];
spf = dir([p2f,'Spike_data_CSC_*.ncs_stimlocked.mat']);

%% tuning log file
f.p2logf = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',pID,'/log/Tuning/'];
f.logf = LF{sesh};
p.ncols = 8;

[LogDat] = getNewLogDataTune(f,p);

%% extract image data tuning
buf = LogDat.dat;%  log{1}.LogDat.dat(2:end,1);
logpicident =zeros(1,length(buf));
parfor n = 1:length(buf)
    logpicident(n) = str2double(buf{n});% each img has a label
end;
numpics = unique(logpicident);

picnamev = cell(1,length(numpics));
for n = 1:length(numpics)
    f = find(logpicident == numpics(n),1,'first');
    picnamev{n} = LogDat.dat{f,2};% name of pic events
end;

titles = cell(1,length(numpics));
parfor p = 1:length(picnamev)
    [a,titles{p},b] = fileparts(picnamev{p});%get the filename of each event
end;

%% settings for kernel smoothing
wsz = 21;
halfwin = floor(wsz/2);
wst = 3;
gw = normpdf(linspace(-halfwin,halfwin,wsz),0,wst);
gw = gw./sum(gw);

%% loop over channels
for it = [30 ]%1:length(spf)%41:48
    
    id1 = regexp(spf(it).name,'CSC')+4;
    id2 = regexp(spf(it).name,'.ncs')-1;
    
    [chan] = spf(it).name(id1:id2);% extract the channel number
    
    muaf = dir([p2f,'MUA_data_CSC_',chan,'.ncs_stimlocked.mat']);
    lfpf = dir([p2f,'LFP_data_CSC_',chan,'.ncs_stimlocked.mat']);
    
    %make data structure with spike times, power, erp, mua, etc
    dat = load([p2f,spf(it).name]);  % spike times  
    dat.mua = load([p2f,muaf.name]); % MUA
    dat.lfp = load([p2f,lfpf.name]); % lfp data (power & erp)
    
    %% erp data processing
    cfg = [];
    cfg.hilbert = 'angle';
    
    [phi] = ft_preprocessing(cfg,dat.lfp.tlck);% extract single trial phase

    n = size(phi.trial,1);
    [itc] = abs(sum(phi.trial,1)./sqrt(n*sum(abs(phi.trial).^2,1)));% compute ITC

    idx = find(itc == max(itc));
    erfm = max(abs(squeeze(dat.lfp.tlck.trial(:,1,idx))),3);% extract ERP at max ITC value
    
    cfg = [];
    cfg.latency = [.1 .9];
    
    [erf] = ft_selectdata(cfg, dat.lfp.tlck );
    
    erf.trial = (erf.trial).^2;%square the amplitude
    
    cfg = [];
    cfg.avgovertime = 'yes';
    
    [erf] = ft_selectdata(cfg, erf );% computes the average of ERP over time
    
    % computes the RMS
    erf.avg = sqrt(erf.avg);
    erf.trial = sqrt(erf.trial);
    %% extract high-gamma power
    cfg = [];
    cfg.latency = [.1 .9];% stim period
    cfg.frequency = [100 200];
    cfg.avgoverfreq = 'yes';
    
    [pow] = ft_selectdata(cfg,dat.lfp.powH);
    
    cfg = [];
    cfg.latency = [-.35 -.075];%baseline period
    cfg.frequency = [100 200];
    cfg.avgoverfreq = 'yes';
    
    [base] = ft_selectdata(cfg,dat.lfp.powH);
    
    b = squeeze(median(mean(base.powspctrm,4),1));
    
    pow.powspctrm = 100*((pow.powspctrm)./b);% remove estimate of baseline activity from single trial HFG 
    
    cfg = [];    
    cfg.avgovertime = 'yes';
    
    [pow] = ft_selectdata(cfg, pow );%computes the avg-HFG over time 
    
    %% extract the MUA
    cfg = [];
    cfg.latency = [.1 .9];%stim period
    
    [mua] = ft_selectdata(cfg,dat.mua.mua);
    
    cfg = [];
    cfg.latency = [-.35 -.075];%baseline period
    
    [base] = ft_selectdata(cfg,dat.mua.mua);
    
    b = squeeze(median(mean(base.trial,3),1));
    
    mua.trial = mua.trial - b;% adjust MUA by substracting baseline estimate from each single trial
        
    cfg = [];    
    cfg.avgovertime = 'yes';
    
    [mua] = ft_selectdata(cfg, mua );% avg-MUA over time
    
    %% extract the spike rate
    cfg                     = [];
    cfg.outputunit          = 'spikecount';
    cfg.keeptrials          = 'yes';
    cfg.binsize             = 0.0006;
    cfg.latency          = [-.5 2.3];% poststim period
    
    [rate]              = ft_spike_psth( cfg,dat.spike3 );%make a psth for each ind trial
    
    for z = 1:size(rate.trial,1)
        rate.trial(z,:,:) = conv(squeeze(rate.trial(z,:,:)),gw,'same');% smooth psth with kernel
        rate.trial(z,:,:) = rate.trial(z,:,:).*rate.fsample;% convert histogram to firing rate
    end;
    
    cfg                     = [];
    cfg.latency             = [-.35 -.075];%baseline period
    cfg.avgovertime         = 'no';
    
    [rateb]                   = ft_selectdata(cfg,rate);% estimate of baseline spiking activity
    
    rate.trial = rate.trial - median(mean(rateb.trial,3),1);% adjust firing rate by substracting baseline activity
    
    cfg                     = [];
    cfg.latency             = [.1 .9];
    cfg.avgovertime         = 'yes';
    
    [rate]                   = ft_selectdata(cfg,rate);% compute average firing rate during stimulation period
    
    %% place all the features into feature-matrix X
    X = zeros(length(numpics),5,2);
    for n = 1:length(numpics)        
        
        f = find(logpicident == numpics(n));% get the trials that correspond to same stim event        
        
        % mean & SD HFG
        X(n,1,1) = squeeze(mean(mean(pow.powspctrm(f,:,:,:),4)));
        X(n,1,2) = squeeze(std(mean(pow.powspctrm(f,:,:,:),4)))/sqrt(length(f));
        % mean & SD MUA
        X(n,2,1) = mean(mua.trial(f));
        X(n,2,2) = std(mua.trial(f))/sqrt(length(f));
        % mean & SD firing rate
        X(n,3,1) = mean(rate.trial(f));
        X(n,3,2) = std(rate.trial(f))/sqrt(length(f));
        % mean & SD MAX of ERP
        X(n,4,1) = mean(erfm(f));
        X(n,4,2) = std(erfm(f))/sqrt(length(f));
        % mean & SD RMS of ERP
        X(n,5,1) = mean(erf.trial(f));
        X(n,5,2) = std(erf.trial(f))/sqrt(length(f));
        
    end;
    %%
    self = find(dat.lfp.powH.freq >=100 & dat.lfp.powH.freq <=200);
    
    selt1 = find(dat.lfp.tlck.time >= 0.1 & dat.lfp.tlck.time <= 0.9);
    selt2 = find(dat.lfp.powH.time >= 0.1 & dat.lfp.powH.time <= 0.9);
    selt3 = find(dat.mua.mua.time >= 0.1 & dat.mua.mua.time <= 0.9);
    
    Y1 = zeros(length(numpics),length(selt1));
    Y2 = zeros(length(numpics),length(selt2));
    Y3 = zeros(length(numpics),length(selt3));
    for kt = 1:length(numpics)
        ix = find(logpicident==numpics(kt));%group data according to stim
        Y1(kt,:)= squeeze(mean(dat.lfp.tlck.trial(ix,:,selt1),1));
        Y2(kt,:) = squeeze(mean(mean(dat.lfp.powH.powspctrm(ix,:,self,selt2),3),1));
        Y3(kt,:) = squeeze(mean(dat.mua.mua.trial(ix,:,selt3),1));
    end;
    
    figure;
    subplot(3,3,1);
    a = gca;
    imagesc(dat.lfp.tlck.time(selt1),1:size(Y1,1),Y1);
    %caxis([min(min(Y))/100*5 max(max(Y))/100*50]);
    title({['channel ',chan];'ERP'});
    
    subplot(3,3,4);
    a = [a gca];
    imagesc(dat.lfp.powH.time(selt2),1:size(Y2,1),Y2);
    %caxis([min(min(Y))/100*5 max(max(Y))/100*50]);
    title(['HFG']);
    
    subplot(3,3,7);
    a = [a gca];
    imagesc(dat.mua.mua.time(selt3),1:size(Y3,1),Y3);
    %caxis([min(min(Y))/100*5 max(max(Y))/100*50]);
    title(['MUA']);
    
    axis(a,'xy');
    %set(a,'XLim',[-.3 1]);
    for kt = 1:length(a)
        hold(a(kt),'on');
        plot(a(kt),[0 0],[1 length(ix)],'w');
        xlabel(a(kt),'Time (s)');
        ylabel(a(kt),'Stim. event (#)');
    end;
    
    subplot(3,3,2:3);
    a = gca;
    errorbar(1:size(X,1),squeeze(X(:,5,1)),squeeze(X(:,5,2)),'ks','MarkerFaceColor','c');
    ylabel('RMS (\muV)');
    
    subplot(3,3,5:6);
    a = [a gca];
    errorbar(1:size(X,1),squeeze(X(:,1,1)),squeeze(X(:,1,2)),'ks','MarkerFaceColor','c');
    ylabel('Relative change (%)');
    
    subplot(3,3,8:9);
    a = [a gca];
    errorbar(1:size(X,1),squeeze(X(:,2,1)),squeeze(X(:,2,2)),'ks','MarkerFaceColor','c');
    ylabel('Relative change (%)');
    
    axis(a,'tight');
        for kt = 1:length(a)
        xlabel(a(kt),'Stim. event (#)');
    end;
    set(gcf,'Color','w');
    %%    
    x6 = [];
    for kt = 1:length(numpics)
        cfg = [];
        cfg.latency = [-.5 1];
        cfg.binsize = 0.051;
        cfg.outputunit = 'rate';
        
        cfg.trials = find(logpicident== numpics(kt));
        psth = ft_spike_psth(cfg,dat.spike3);
        
        x6(kt,:) = psth.avg;
    end;
    
    figure;
    subplot(4,1,1:3);
    imagesc(psth.time,1:size(x6,1),x6);
    xlim([.1 .9]);axis xy;
    title({['channel ',chan];'PSTH'});
    xlabel('Time (s)');
    ylabel('Stim event (#)');
    
    subplot(4,1,4);
    errorbar(1:size(X,1),squeeze(X(:,3,1)),squeeze(X(:,3,2)),'ks','MarkerFaceColor','k');
    axis tight;
    xlabel('Stim event (#)');
    ylabel('Spikes/s');
    set(gcf,'Color','w');
    %%
    [r,p] = corr(squeeze(X(:,:,1)),'Type','Spearman');
    figure;
    imagesc(r);
    set(gca,'XTick',1:5);
    set(gca,'YTick',1:5);
    set(gca,'YTickLabel',{'HFG','MUA','FR','MAX','RMS'});
    set(gca,'XTickLabel',{'HFG','MUA','FR','MAX','RMS'});
    caxis([-1 1]);
    colormap jet;
    title({['channel ',chan];'Correlation matrix'});
end;