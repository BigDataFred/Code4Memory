%%
restoredefaultpath;
addpath('~rouxf/fieldtrip-20161009/');
addpath(genpath('~rouxf/AnalysisFred/'));

ft_defaults;
%%
if isempty(gcp('nocreate'))
    parpool(32);
end;
%%
p2f = '/home/rouxf/res/tuning/P02/2016-07-09_11-13-28/';
files = dir([p2f,'Spike_data_CSC_*_stimlocked.mat']);
%%
f = [];
f.p2logf = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/P02/log/Tuning/'];
f.logf = 'P02_TS01_log_ctune_09-Jul-2016_11_40_xx.txt';
p.ncols = 8;

[LogDat] = getNewLogDataTune(f,p);
%%
ch = cell(1,length(files)); 

parfor it = 1:length(files) 
    
    id1 = regexp(files(it).name,'CSC')+5;
    id2 = regexp(files(it).name,'.ncs')-1;
    
    ch{it} =  str2double( files(it).name(id1:id2) );
    if isnan(ch{it})
        id1 = regexp(files(it).name,'CSC')+4;
        id2 = regexp(files(it).name,'.ncs')-1;
        ch{it} = files(it).name(id1:id2);
    end;
end;
%%

[ch,s_idx] = sort(ch);

files  = files(s_idx);
%%
dat = cell(1,length(files));
parfor it = 1:length(files)
    
    dat{it} = load([p2f,files(it).name]);
end;
%%
chanID = ( ch );%makeCSClabels();
n = length(chanID)/8;
%%
buf = LogDat.dat;%  log{1}.LogDat.dat(2:end,1);
logpicident =[];
parfor n = 1:length(buf)
    logpicident(n) = str2num(buf{n});
end;
numpics = length(unique(logpicident));

picnamev = [];
for n = 1:numpics
    f = find(logpicident == n,1,'first');
    picnamev{n} = LogDat.dat{f,2};
end;

titles = {};
parfor p = 1:length(picnamev)
    [a,titles{p},b] = fileparts(picnamev{p});
end;
%%
psth = cell(1,length(dat));
x = zeros(length(dat),115);
parfor it = 1:length(dat)
    if isfield(dat{it}.spike3,'cfg')
        dat{it}.spike3 = rmfield(dat{it}.spike3,'cfg');
    end;
    
    cfg                     = [];
    cfg.binsize             = 0.01;
    cfg.outputunit          = 'rate';
    cfg.latency             = [-0.15 1];
    
    [psth{it}]              = ft_spike_psth(cfg, dat{it}.spike2);
    if isfield(dat{it}.spike3,'cfg')
        dat{it}.spike3 = rmfield(dat{it}.spike2,'cfg');
    end;
    
    dum = psth{it}.avg;
    
    b = dum(psth{it}.time <0);
    sd = repmat(std(b,0,2),[1 size(dum,2)]);
    dum(find(psth{it}.time >= -0.005 & psth{it}.time <= 0.005)) = NaN;
    dum = dum - repmat(mean(b),[1 length(dum)]);
    
    x(it,:) = dum;
    
    [dat{it}.spike3] = sortST2RT(LogDat.RT,dat{it}.spike3);
    
end;
z = max(x,[],2);
z = (z-mean(z))./std(z);

sel_idx = find(z>2);
z2 = (x-min(min(x)))./(max(max(x)) - min(min(x)));
%%
ntrl = length(unique(dat{1}.spike3.trial{1}));
n = length(chanID)/8;

figure;
for it = 1:length(dat)
    
    
    subplot(n,8,it);

    hold on;
    
    chck = (dat{it}.spike3.time{1} >= -.15 & dat{it}.spike3.time{1} <=1);
    chck2 = ismember(dat{it}.spike3.trial{1},1:ntrl);
    x = dat{it}.spike3.time{1}(chck & chck2);
    y = dat{it}.spike3.trial2{1}(chck & chck2);
    x = x(:)';
    y = y(:)';
    y = [y-0.45;y+0.45];
    x = [x;x];
    
    line(x,y,'Color','k');
    %plot(x(1,:),y(1,:),'k.');
    
    plot([0 0],[1 ntrl],'r');
        
    plot(psth{it}.time,z2(it,:).*ntrl,'b');
    
    if ismember(it,sel_idx)
        box(gca,'on');
        set(gca,'LineWidth',3);
        set(get(gca,'Xaxis'),'Color',[0 0 .9]);
        set(get(gca,'Yaxis'),'Color',[0 0 .9]);
    end;
    axis tight;
    xlim([-.15 1]);
    set(gca,'YDir','normal');
    xlabel('Time (s)');
    ylabel('Trial #');
    title(chanID(it));
end;
set(gcf,'Color','w');
%%
ntrl = length(unique(dat{1}.spike3.trial{1}));

for it = 1:length(sel_idx)
    figure;
    
    hold on;
    
    chck = (dat{sel_idx(it)}.spike3.time{1} >= -.15 & dat{sel_idx(it)}.spike3.time{1} <=1);
    chck2 = ismember(dat{sel_idx(it)}.spike3.trial{1},1:ntrl);
    x = dat{sel_idx(it)}.spike3.time{1}(chck & chck2);
    y = dat{sel_idx(it)}.spike3.trial2{1}(chck & chck2);
    x = x(:)';
    y = y(:)';
    y = [y-0.45;y+0.45];
    x = [x;x];
    
    line(x,y,'Color','k');
    
    plot(sort(LogDat.RT)-1,1:length(LogDat.RT),'r.');
    
    plot([0 0],[1 ntrl],'r');
    
    plot(psth{sel_idx(it)}.time,z2(sel_idx(it),:).*ntrl,'b');
    
    box(gca,'on');
    set(gca,'LineWidth',3);
    set(get(gca,'Xaxis'),'Color',[0 0 .9]);
    set(get(gca,'Yaxis'),'Color',[0 0 .9]);
    axis tight;
    xlim([-.15 1]);
    set(gca,'YDir','normal');
    xlabel('Time (s)');
    ylabel('Trial #');
    title(chanID(sel_idx(it)));
end;
set(gcf,'Color','w');

%%
wsz = 21;
halfwin = floor(wsz/2);
wst = 3;
gw = normpdf(linspace(-halfwin,halfwin,wsz),0,wst);
gw = gw./sum(gw);
%%
r = cell(1,length(dat));
s = cell(1,length(dat));

for jt = 1:length(dat)

        
        r{jt} = zeros(length(LogDat.ID),2);
        s{jt} = zeros(length(LogDat.ID),1);
        po = [];
 
        cfg                     = [];
        cfg.outputunit          = 'spikecount';
        cfg.keeptrials          = 'yes';
        cfg.binsize             = 0.0006;
        cfg.latency          = [-.5 2.3];
        
        [rate]              = ft_spike_psth( cfg,dat{jt}.spike2 );
        
        for z = 1:size(rate.trial,1)
            rate.trial(z,:,:) = conv(squeeze(rate.trial(z,:,:)),gw,'same');
            rate.trial(z,:,:) = rate.trial(z,:,:).*rate.fsample;
        end;
        
        cfg                     = [];
        cfg.latency             = [-.35 -.075];
        cfg.avgovertime         = 'no';
        
        [rateb]                   = ft_selectdata(cfg,rate);
        
        rate.trial = rate.trial - median(mean(rateb.trial,3),1);
        
        cfg                     = [];
        cfg.latency             = [.1 .9];
        cfg.avgovertime         = 'yes';
        
        [rate]                   = ft_selectdata(cfg,rate);
        %%
        for p = 1:numpics
                      
            f = find(logpicident == p);
                                    
            r{jt}(p,:) = [mean(rate.trial(f)) std(rate.trial(f))/sqrt(length(f))];
                        
            [h,pval] = ttest(rate.trial(f),0,'tail','right');
            
            if h==1 && pval < 0.001%
                s{jt}(p) = 1;
            end;
            
        end;    
end;
%%
n = length(ch)/8;

figure;
lm = zeros(length(r),2);
a = zeros(length(r),1);
xl = cell(1,length(r));
for jt = 1:length(r)
    
    sel_idx = [];    
        
    subplot(n,8,jt);
    hold on;
    a(jt) = gca;
    errorbar(1:size(r{jt},1),r{jt}(:,1),r{jt}(:,2),'s','Color',[.75 .75 .75]);
    
    chck = s{jt};
    sel_idx = find(chck ==1);
    
    if ~isempty(sel_idx)
        errorbar(sel_idx,r{jt}(sel_idx,1),r{jt}(sel_idx,2),'ks');
        plot(sel_idx,r{jt}(sel_idx,1),'rs','MarkerFaCeColor','r');
    end;
    
    lm(jt,:) = [min(r{jt}(:,1)-r{jt}(:,2)) max(r{jt}(:,1)+r{jt}(:,2))];
    
    xlabel('Stim.');
    ylabel('(spikes/sec)');
    title(chanID(jt));
    set(gca,'XTick',sel_idx);
    %set(gca,'XTickLabel',titles(sel_idx));
    set(gca,'XTickLabelRotation',65);
    set(gca,'YTick',[min(r{jt}(:,1)-r{jt}(:,2)) max(r{jt}(:,1)+r{jt}(:,2))]);
    set(gca,'YTickLabel',round(get(gca,'YTick')*10)/10);
    sel_idx = [];
    axis(gca,'tight');
end;

 
delete(gcp);