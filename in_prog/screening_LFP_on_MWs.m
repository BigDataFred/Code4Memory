%%
addpath('/media/rouxf/rds-share/Common/fieldtrip-20170115/');
ft_defaults;
addpath('/media/rouxf/rds-share/Fred/code/mcode/custom/project_EM/');

%%
[rpath] = '/media/rouxf/rds-share/iEEG_DATA/MICRO/P02/fVSpEM/';

%% check for existing sessions
chck = dir(rpath);
chck(1:2) = [];
ix = [];
c =0;
for it = 1:length(chck)
    if chck(it).isdir
        c =c+1;
        ix(c) = it;
    end;
end;

sesh = chck(ix);

if isempty(sesh)
    error('there is no corresponding data');
end;

%% loop over sessions and extract behavioral data

for nt = 1:length(sesh)
    
    p2log = [rpath,sesh(nt).name,filesep,'log_dat',filesep];
    log_fn = dir([p2log,'*LogDat.mat']);
    
    load([p2log,log_fn.name]);% load the data
    
    p2lfp = [rpath,sesh(nt).name,filesep,'lfp_dat',filesep];
    lfp_fn = dir([p2lfp,'*lfp_dat*.mat']);
    
    dum = cell(1,length(lfp_fn));
    parfor it = 1:length(lfp_fn)
        fprintf([num2str(it),'/',num2str(length(lfp_fn))]);
        dat = load([p2lfp,lfp_fn(it).name]);
        dum{it} = dat.save_data{1}{1}{1};
        fprintf('\n');
    end;
    
    [dum] = ft_appenddata([],dum{:});
    
    [trl_all] = extract_Enc_trl(LogDat1);% all, remembered (r), not-remembered (nr)
    
    lfp_dat = dum;
    
    %%
    cfg                     = [];
    cfg.bpfilter            = 'yes';
    cfg.bpfreq              = [0.5 30];
    
    [dum] = ft_preprocessing( cfg , dum );
    
    cfg                     = [];
    cfg.trials              = trl_all;
    cfg.latency              = [-2 5];
    
    [dum] = ft_selectdata( cfg , dum );
    
    %%
    cfg                     = [];
    cfg.keeptrials          = 'yes';
    
    [erf] = ft_timelockanalysis( cfg , dum );
    
    cfg                     = [];
    cfg.baseline            = [-1 -.1];
    
    [erf] = ft_timelockbaseline( cfg, erf );
    
    %%
    cfg                     = [];
    cfg.hilbert             = 'complex';
    
    [dum2] = ft_preprocessing( cfg , dum );
    
    cfg                     =  [];
    cfg.keeptrials          = 'yes';
    
    [dum2] = ft_timelockanalysis( cfg , dum2 );
    
    [itc] = dum2.trial./abs(dum2.trial);
    N  = size(itc,1);
    itc = sum(itc,1);
    itc =  abs( itc )./N;
    itc = squeeze(itc);
    
    %%
    nChans = length(dum.label)/8;
    
    selIx = 0:8:(8*nChans);
    
    for it = 1:length( selIx)-1
        
        ix = selIx(it)+1:selIx(it+1);
        
        chanLabel = erf.label{ix(1)};
        chanLabel(regexp(chanLabel,'\d{1}')) = [];
        chanLabel(regexp(chanLabel,'_')) = [];
        
        cfg                     = [];
        cfg.channel             = ix;
        cfg.avgoverchan         = 'no';
        
        [avg] = ft_selectdata( cfg , erf );
        
        figure;
        subplot(4,1,1:3);
        hold on;
        imagesc(avg.time,1:size(avg.trial,1),squeeze(mean(avg.trial,2)));
        plot([0 0],[1 size(avg.trial,1)],'w-');
        plot([2 2],[1 size(avg.trial,1)],'w-');
        axis xy;axis tight;
        title(chanLabel);
        
        subplot(4,1,4);
        hold on;
        Y= avg.avg;
        plot(avg.time,Y,'k');
        plot([0 0],[min(min(avg.avg)) max(max(avg.avg))],'r-');
        plot([2 2],[min(min(avg.avg)) max(max(avg.avg))],'r-');
        axis tight;
        
        
    end;
    
end;

% %%
% figure;
% for it = 1:length(erf.label)
%     six = it;
%
%     subplot(311);
%     hold on;
%     imagesc(erf.time,1:size(erf.trial,1),squeeze(erf.trial(:,six(end),:)));
%     plot([0 0],[1 size(erf.trial,1)],'w');axis xy;axis tight;
%     plot([2 2],[1 size(erf.trial,1)],'w');axis xy;axis tight;
%     %set(gca,'Xlim',[1.96 2]);
%
%     subplot(312);
%     hold on;
%     Y = squeeze(erf.trial(:,six(end),:));
%     plot(erf.time,squeeze(Y),'k');
%     plot(erf.time,mean(Y,1),'Color','g');
%     plot([0 0],[min(min(Y)) max(max(Y))],'r');
%     plot([2 2],[min(min(Y)) max(max(Y))],'r');
%     axis(gca,'tight');
%     %set(gca,'Xlim',[1.96 2]);
%
%     subplot(313);
%     hold on;
%     Y = itc(six(end),:);
%     plot(erf.time,Y,'b');
%     plot([0 0],[0 1],'r');
%     plot([2 2],[0 1],'r');
%     axis tight;
%     %set(gca,'Xlim',[1.96 2]);set(gca,'Ylim',[0 1]);
%
%     pause;
%     clf;
%
% %     %%
% %     cfg                     = [];
% %     cfg.detrend             = 'yes';
% %     cfg.demean              = 'yes';
% %
% %     [lfp_dat] = ft_preprocessing( cfg , lfp_dat );
% %
% %     %%
% %     cfg                     = [];
% %     cfg.method              = 'mtmconvol';
% %     cfg.pad                 = 'maxperlen';
% %     cfg.taper               = 'dpss';
% %     cfg.tapsmofrq           = 1;
% %     cfg.foi                 = 1:25;
% %     cfg.t_ftimwin           = ones(1,length(cfg.foi));
% %     cfg.keeptrials          = 'yes';
% %     cfg.toi                 = lfp_dat.time{1}(1):0.065:lfp_dat.time{1}(end);
% %
% %     [powL] = ft_freqanalysis( cfg , lfp_dat );
% %
% %
% %     cfg                     = [];
% %     cfg.method              = 'mtmconvol';
% %     cfg.pad                 = 'maxperlen';
% %     cfg.taper               = 'dpss';
% %     cfg.tapsmofrq           = 10;
% %     cfg.foi                 = 25:170;
% %     cfg.t_ftimwin           = 0.5*ones(1,length(cfg.foi));
% %     cfg.keeptrials          = 'yes';
% %     cfg.toi                 = lfp_dat.time{1}(1):0.065:lfp_dat.time{1}(end);
% %
% %     [powH] = ft_freqanalysis( cfg , lfp_dat );
% %
% %     %%
% %     Y1 = squeeze(mean(powL.powspctrm(ix{4},:,:,:),1));
% %     Y2 = squeeze(mean(powL.powspctrm([ix{5};ix{6}],:,:,:),1));
% %
% %     figure;
% %     imagesc(powL.time,powL.freq,squeeze(mean(Y2(2:end,:,:)-Y1(2:end,:,:),1)));
% %     axis xy;
%
% end;

%%

