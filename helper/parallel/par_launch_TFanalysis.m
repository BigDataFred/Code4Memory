%%
restoredefaultpath;
addpath('~rouxf/fieldtrip-20161009/');
ft_defaults;
addpath(genpath('/media/rouxf/rds-share/iEEG_DATA/MICRO/analysis_pipeline/'));

%%
p2d = '/media/rouxf/rds-share/iEEG_DATA/MICRO/P04/fvsPEM/';
d = dir([p2d,'2016-*']);

%%
y = [];
for it = 1:length(d)
    f = dir([p2d,d(it).name , filesep , 'log_dat', filesep , '*_EMtask_LogDat.mat']);
    load([ p2d , d(it).name , filesep , 'log_dat', filesep , f.name ],'ix','ix_readme','RTs');
    y = [y;RTs];
end;

%%
hits = {};  misses = {};    
base = {};  post = {};
for it = 1:length(d)
    
    f = dir([p2d,d(it).name , filesep , 'lfp_dat', filesep , '*_EMtask_preprocLFPdata_ds.mat']);
    load([ p2d , d(it).name , filesep , 'lfp_dat', filesep , f.name ]);
    
    f = dir([p2d,d(it).name , filesep , 'log_dat', filesep , '*_EMtask_LogDat.mat']);
    load([ p2d , d(it).name , filesep , 'log_dat', filesep , f.name ],'ix','ix_readme','RTs');   
        
    ch = {};
    for zt = 1:length(lfp_data.label)
        ch{zt} = lfp_data.label{zt}(regexp(lfp_data.label{zt},'_')+1:end);
    end;
    
    chck1 = regexp(ch,'L');
    chck2 = regexp(ch,'R');
    
    c = 0;sel_idx1 = [];
    for zt = 1:length(chck1)
        if ~isempty(chck1{zt})
            c = c+1;
            sel_idx1(c) = zt;
        end;
    end;
    
        c = 0;sel_idx2 = [];
    for zt = 1:length(chck2)
        if ~isempty(chck2{zt})
            c = c+1;
            sel_idx2(c) = zt;
        end;
    end;
    
    ch2 = [];
    for zt = 1:length(ch)
        ch2(zt) = str2double(ch{zt}(2:end));
    end;
    
    [~,s_idx1] = sort(ch2(sel_idx1));
    [~,s_idx2] = sort(ch2(sel_idx2));
    
    sel_idx1 = sel_idx1(s_idx1);
    sel_idx2 = sel_idx2(s_idx2);
    
    s_idx = [sel_idx1 sel_idx2]';
    
    
    lfp_data.label = lfp_data.label(s_idx);
    for zt = 1:length(lfp_data.trial)
        lfp_data.trial{zt} = lfp_data.trial{zt}(s_idx,:);
    end;
    
    cfg                     = [];
    cfg.latency             = [2 3+floor(min(y))];
    cfg.channel             = lfp_data.label(1:24);
    
    [sel] = ft_selectdata( cfg, lfp_data );
    
    n = length(sel.time{1});
    
    [pow] = compute_LFP_spectrum( 'low', sel , 8 );
    
    post{it} = pow{1};
    
    cfg                     = [];
    cfg.trials              = [ix{4}];
    
    [hits{it}] = ft_selectdata( cfg, pow{1}  );
    
    cfg                     = [];
    cfg.trials              = [ix{6}];
    
    [misses{it}] = ft_selectdata( cfg, pow{1}  );
    
    cfg                     = [];
    cfg.latency             = [-.75 -.25];
    cfg.channel             = lfp_data.label(1:24);
    
    [sel] = ft_selectdata( cfg, lfp_data );
    
    n2 = length(sel.time{1});
    f = (n-n2)/2;
    
    pad = zeros(length(sel.label),f);
    
    for kt = 1:length(sel.trial)
        sel.trial{kt}   = [pad sel.trial{kt} pad];
        sel.time{kt}    = 0:1/sel.fsample:(length(sel.trial{kt})-1)/sel.fsample;
    end;
    n2 = length(sel.time{1});
    
    if n2 ~= n
        error('baseline and poststim periods must match');
    end;
    
    [pow] = compute_LFP_spectrum( 'low', sel , 8 );
    
    base{it} = pow{1};
end;
%%
cfg                     = [];
cfg.parameter           = 'powspctrm';

[base] = ft_appendfreq( cfg , base{:} );

[post] = ft_appendfreq( cfg , post{:} );

[hits] = ft_appendfreq( cfg , hits{:} );

[misses] = ft_appendfreq( cfg , misses{:} );

ntrl = size(post.powspctrm,1);
ntrl1 = size(hits.powspctrm,1);
ntrl2 = size(misses.powspctrm,1);

M =  repmat(mean(base.powspctrm,1),[ntrl 1 1]);
SD = repmat(std(base.powspctrm,0,1),[ntrl 1 1]);

M1 =  repmat(mean(base.powspctrm(ix{4},:,:),1),[ntrl1 1 1]);
SD1 = repmat(std(base.powspctrm(ix{4},:,:),0,1),[ntrl1 1 1]);

M2  = repmat(mean(base.powspctrm(ix{6},:,:),1),[ntrl2 1  1]);
SD2 = repmat(std(base.powspctrm(ix{6},:,:),0,1),[ntrl2 1 1]);
%%
foi = [ 6 8 ]
PTA1 = cell(1,2);
PTA2 = cell(1,2);
for it = 1:length(d)
    
    f = dir([p2d,d(it).name , filesep , 'lfp_dat', filesep , '*_EMtask_preprocLFPdata_ds.mat']);
    load([ p2d , d(it).name , filesep , 'lfp_dat', filesep , f.name ]);
          
    cfg                     = [];
    cfg.latency             = [-.75 -.25];
    cfg.channel = lfp_data.label([11 13]);
    
    [sel] = ft_selectdata( cfg, lfp_data );
    
    cfg = [];    
    [sel2] = ft_selectdata(cfg,sel);
    
    cfg = [];    
    cfg.bpfilter = 'yes';
    cfg.bpfreq = foi;
    cfg.bfilttype = 'but';
    cfg.padding = 10;
    cfg.padtype = 'mirror';
    
    [bpf] = ft_preprocessing(cfg,sel);
    
    cfg = [];
    cfg.keeptrials = 'yes';
    
    [bpf] = ft_timelockanalysis(cfg,bpf);
    [sel2] = ft_timelockanalysis(cfg,sel2);
    
    for ft = 1:length(sel2.label)
    
    c = 0;pTA1 = [];
    ntrl =size(bpf.trial,1);    
    for kt = 1:ntrl
    sig = squeeze(bpf.trial(kt,ft,:));
    k =  local_max(sig);
    k(find(diff(k) <(floor(1/foi(2)*sel.fsample))))=[];
    
    sig = squeeze(sel2.trial(kt,ft,:));
    win = 2*round(1/foi(1)*sel.fsample/2); 
    
    for gt = 1:length(k)
        
        if k(gt)-win >0 && k(gt)+win < length(sig)
            c = c+1;
            pTA1(c,:) =sig(k(gt)-win:k(gt)+win)';
        end;
        
    end;    
    end;
    PTA1{ft} = [PTA1{ft};pTA1];    
    end;
    
        cfg                     = [];
    cfg.latency             = [2 3+floor(min(y))];
    cfg.channel = lfp_data.label([11 13]);
    
    [sel] = ft_selectdata( cfg, lfp_data );
    
    cfg = [];    
    [sel2] = ft_selectdata(cfg,sel);
    
    cfg = [];
    cfg.bpfilter = 'yes';
    cfg.bpfreq = foi;
    cfg.bfilttype = 'but';
    cfg.padding = 10;
    cfg.padtype = 'mirror';
    
    [bpf] = ft_preprocessing(cfg,sel);
    
    cfg = [];
    cfg.keeptrials = 'yes';
    
    [bpf] = ft_timelockanalysis(cfg,bpf);
    [sel2] = ft_timelockanalysis(cfg,sel2);
    
    pTA2 = [];
    c = 0;
    ntrl =size(bpf.trial,1);
    for ft= 1:length(sel2.label)
    for kt = 1:ntrl
    sig = squeeze(bpf.trial(kt,ft,:));
    k =  local_max(sig);
    k(find(diff(k) <(floor(1/foi(2)*sel.fsample))))=[];
    
    sig = squeeze(sel2.trial(kt,ft,:));
    win = 2*round(1/foi(1)*sel.fsample/2); 
    
    for gt = 1:length(k)
        
        if k(gt)-win >0 && k(gt)+win < length(sig)
            c = c+1;
            pTA2(c,:) =sig(k(gt)-win:k(gt)+win)';
        end;
        
    end;
    end;
    
    PTA2{ft} = [PTA2{ft};pTA2];
    end;
    
end;

dt = win/sel.fsample;
n = win*2+1;
t = linspace(-dt,dt,n);

figure;
subplot(211);
hold on;
plot(t,mean(PTA1{1},1));
plot(t,mean(PTA2{1},1),'r');
%n = size(PTA1{1},1);
%errorbar(t,mean(PTA1{1},1),std(PTA1{1},1)./sqrt(n-1),'bs');
%n = size(PTA2{1},1);
%errorbar(t,mean(PTA2{1},1),std(PTA2{1},1)./sqrt(n-1),'rs');
axis tight;

subplot(212);
hold on;
plot(t,mean(PTA1{2},1));
plot(t,mean(PTA2{2},1),'r');
axis tight;
%%
dum = [];
dum = [squeeze(mean(base.powspctrm,1));squeeze(mean(post.powspctrm,1))];

lm = [];
lm(1) = min(min( dum ));
lm(2) = max(max( dum ));

ntrl = size(post.powspctrm,1);
figure;
for it = 1:24
    
M1 = squeeze(mean((base.powspctrm(:,it,:)),1));
M2 = squeeze(mean((post.powspctrm(:,it,:)),1));

SD1 = squeeze(std((base.powspctrm(:,it,:)),1))./sqrt(ntrl-1);
SD2 = squeeze(std((post.powspctrm(:,it,:)),1))./sqrt(ntrl1);

    subplot(24/8,8,it);
    hold on;
    plot(base.freq,M1,'Color',[.75 .25 .25],'LineWidth',2);
    plot(base.freq,M1-SD1,'Color',[.75 .25 .25],'LineWidth',1);
    plot(base.freq,M1+SD1,'Color',[.75 .25 .25],'LineWidth',1);
    
    plot(post.freq,M2,'Color',[0 0 .5],'LineWidth',2);
    plot(post.freq,M2-SD2,'Color',[0 0 .5],'LineWidth',2);
    plot(post.freq,M2+SD2,'Color',[0 0 .5],'LineWidth',2);
    
    axis tight;
    ylim([lm]);
    lm2= get(gca,'YLim');
    plot([7.5 7.5],[lm2(1) lm2(2)],'Color',[.25 .25 .25]);
    
    set(gca,'XTick',7.5.*[1:4]);
end;

set(gcf,'Color','w');
%%
idx = 1:49;
M1 = squeeze(mean(post.powspctrm(idx,19,:),1));
SD1 = squeeze(std(post.powspctrm(idx,19,:),1))/sqrt(length(idx)-1);

idx = 50:99;
M2 = squeeze(mean(post.powspctrm(idx,19,:),1));
SD2 = squeeze(std(post.powspctrm(idx,19,:),1))/sqrt(length(idx)-1);

idx = 100:147;
M3 = squeeze(mean(post.powspctrm(idx,19,:),1));
SD3 = squeeze(std(post.powspctrm(idx,19,:),1))/sqrt(length(idx)-1);

idx = 1:49;
M4 = squeeze(mean(post.powspctrm(idx,20,:),1));
SD4 = squeeze(std(post.powspctrm(idx,20,:),1))/sqrt(length(idx)-1);

idx = 50:99;
M5 = squeeze(mean(post.powspctrm(idx,20,:),1));
SD5 = squeeze(std(post.powspctrm(idx,20,:),1))/sqrt(length(idx)-1);

idx = 100:147;
M6 = squeeze(mean(post.powspctrm(idx,20,:),1));
SD6 = squeeze(std(post.powspctrm(idx,20,:),1))/sqrt(length(idx)-1);


figure;
subplot(121);
hold on;
errorbar(post.freq,M1,SD1);
errorbar(post.freq,M2,SD2,'r');
errorbar(post.freq,M3,SD3,'k');
axis tight;
subplot(122);
hold on;
errorbar(post.freq,M4,SD4);
errorbar(post.freq,M5,SD5,'r');
errorbar(post.freq,M6,SD6,'k');
axis tight;
%%
cfg                     = [];
cfg.channel             = post.label([19 20]);
cfg.frequency           = [30 60];%[13 17];%
cfg.avgoverfreq         = 'yes';

[pow1] = ft_selectdata( cfg, post );

cfg                     = [];
cfg.channel             = base.label([19 20]);
cfg.frequency           = [30 60];%[13 17];%
cfg.avgoverfreq         = 'yes';

[pow2] = ft_selectdata( cfg, base );

figure;
subplot(221);
Mx = [];
Mx(1) = mean(pow1.powspctrm(:,1));
Mx(2) = mean(pow1.powspctrm(:,2));
Mx(3) = mean(pow2.powspctrm(:,1));
Mx(4) = mean(pow2.powspctrm(:,2));

SDx = [];
SDx(1) = std(pow1.powspctrm(:,1))/sqrt(ntrl1-1);
SDx(2) = std(pow1.powspctrm(:,2))/sqrt(ntrl1-1);
SDx(3) = std(pow2.powspctrm(:,1))/sqrt(ntrl2-1);
SDx(4) = std(pow2.powspctrm(:,2))/sqrt(ntrl2-1);

hold on;
bar(1,Mx(1),'FaceColor', [.75 .25 .25]);
bar(2,Mx(3),'FaceColor',[0 0 .5]);
bar(3,Mx(2),'FaceColor',[.75 .25 .25]);
bar(4,Mx(4),'FaceColor',[0 0 .5]);

plot([1 1],[Mx(1)-SDx(1) Mx(1)+SDx(1)],'k','LineWidth',3);
plot([2 2],[Mx(3)-SDx(3) Mx(3)+SDx(3)],'k','LineWidth',3);
plot([3 3],[Mx(2)-SDx(2) Mx(2)+SDx(2)],'k','LineWidth',3);
plot([4 4],[Mx(4)-SDx(4) Mx(4)+SDx(4)],'k','LineWidth',3);
plot([.9 1.1],[Mx(1)-SDx(1) Mx(1)-SDx(1)],'k','LineWidth',3);
plot([.9 1.1],[Mx(1)+SDx(1) Mx(1)+SDx(1)],'k','LineWidth',3);
plot([1.9 2.1],[Mx(3)-SDx(3) Mx(3)-SDx(3)],'k','LineWidth',3);
plot([1.9 2.1],[Mx(3)+SDx(3) Mx(3)+SDx(3)],'k','LineWidth',3);
plot([2.9 3.1],[Mx(2)-SDx(2) Mx(2)-SDx(2)],'k','LineWidth',3);
plot([2.9 3.1],[Mx(2)+SDx(2) Mx(2)+SDx(2)],'k','LineWidth',3);
plot([3.9 4.1],[Mx(4)-SDx(4) Mx(4)-SDx(4)],'k','LineWidth',3);
plot([3.9 4.1],[Mx(4)+SDx(4) Mx(4)+SDx(4)],'k','LineWidth',3);

subplot(222);
[n1,x1] = hist(log10([pow1.powspctrm(:,1);pow1.powspctrm(:,2)]));
[n2,x2] = hist(log10([pow2.powspctrm(:,1);pow2.powspctrm(:,2)]));

hold on;
plot(x1,n1./sum(n1),'Color',[.75 .25 .25],'LineWidth',3);
plot(x2,n2./sum(n1),'Color',[0 0 .5],'LineWidth',3);
axis tight;

subplot(223);
Mx = [];
Mx(1) = mean(pow1.powspctrm(ix{4},1));
Mx(2) = mean(pow1.powspctrm(ix{4},2));
Mx(3) = mean(pow1.powspctrm(ix{6},1));
Mx(4) = mean(pow1.powspctrm(ix{6},2));

SDx = [];
SDx(1) = std(pow1.powspctrm(ix{4},1))/sqrt(length(ix{4})-1);
SDx(2) = std(pow1.powspctrm(ix{4},2))/sqrt(length(ix{4})-1);
SDx(3) = std(pow1.powspctrm(ix{6},1))/sqrt(length(ix{6})-1);
SDx(4) = std(pow1.powspctrm(ix{6},2))/sqrt(length(ix{6})-1);

hold on;
bar(1,Mx(1),'FaceColor', [.75 .25 .25]);
bar(2,Mx(3),'FaceColor',[0 0 .5]);
bar(3,Mx(2),'FaceColor',[.75 .25 .25]);
bar(4,Mx(4),'FaceColor',[0 0 .5]);

plot([1 1],[Mx(1)-SDx(1) Mx(1)+SDx(1)],'k','LineWidth',3);
plot([2 2],[Mx(3)-SDx(3) Mx(3)+SDx(3)],'k','LineWidth',3);
plot([3 3],[Mx(2)-SDx(2) Mx(2)+SDx(2)],'k','LineWidth',3);
plot([4 4],[Mx(4)-SDx(4) Mx(4)+SDx(4)],'k','LineWidth',3);

plot([.9 1.1],[Mx(1)-SDx(1) Mx(1)-SDx(1)],'k','LineWidth',3);
plot([.9 1.1],[Mx(1)+SDx(1) Mx(1)+SDx(1)],'k','LineWidth',3);
plot([1.9 2.1],[Mx(3)-SDx(3) Mx(3)-SDx(3)],'k','LineWidth',3);
plot([1.9 2.1],[Mx(3)+SDx(3) Mx(3)+SDx(3)],'k','LineWidth',3);
plot([2.9 3.1],[Mx(2)-SDx(2) Mx(2)-SDx(2)],'k','LineWidth',3);
plot([2.9 3.1],[Mx(2)+SDx(2) Mx(2)+SDx(2)],'k','LineWidth',3);
plot([3.9 4.1],[Mx(4)-SDx(4) Mx(4)-SDx(4)],'k','LineWidth',3);
plot([3.9 4.1],[Mx(4)+SDx(4) Mx(4)+SDx(4)],'k','LineWidth',3);

subplot(224);
[n1,x1] = hist(log10([pow1.powspctrm(ix{4},1);pow1.powspctrm(ix{4},2)]));
[n2,x2] = hist(log10([pow1.powspctrm(ix{6},1);pow1.powspctrm(ix{6},2)]));

hold on;
plot(x1,n1./sum(n1),'Color',[.75 .25 .25],'LineWidth',3);
plot(x2,n2./sum(n1),'Color',[0 0 .5],'LineWidth',3);
axis tight;
%%
dum = [];
dum = squeeze(mean((post.powspctrm-M)./M,1));

lm = [];
lm(1) = min(min( dum ));
lm(2) = max(max( dum ));

figure;
for it = 1:24
    subplot(24/8,8,it);
    hold on;
    plot(base.freq,squeeze(mean((post.powspctrm(:,it,:)-M(:,it,:))./SD(:,it,:),1)),'b');
    axis tight;
    ylim([lm]);
end;
%%
lm = [];
lm(1) = min(min(mean(M1,1)));
lm(2) = max(max(mean(M2,1)));

figure;
for it = 1:24
    subplot(24/8,8,it);
    hold on;
    plot(base.freq,squeeze(mean(M1(:,it,:),1)),'r');
    plot(base.freq,squeeze(mean(M2(:,it,:),1)),'b');

    axis tight;
    ylim([lm]);
end;

%%
lm = [];
lm(1) = min(min(mean(SD1,1)));
lm(2) = max(max(mean(SD2,1)));

figure;
for it = 1:24
    subplot(24/8,8,it);
    hold on;
    plot(base.freq,squeeze(mean(SD1(:,it,:),1)),'r');
    plot(base.freq,squeeze(mean(SD2(:,it,:),1)),'b');

    axis tight;
    ylim([lm]);
end;

%%
dum1 = squeeze(mean(hits.powspctrm,1));
dum2 = squeeze(mean(misses.powspctrm,1));

lm = [];
lm(1) = min(min( [dum1;dum2] ));
lm(2) = max(max( [dum1;dum2] ));

figure;
for it = 1:24
    subplot(24/8,8,it);
    hold on;
    plot(misses.freq,squeeze(mean(misses.powspctrm(:,it,:),1)),'k');
    plot(hits.freq,squeeze(mean(hits.powspctrm(:,it,:),1)),'r');
    axis tight;
    ylim([lm]);
    set(gca,'YTick',[lm(1) lm(2)]);
    set(gca,'YTickLabel',round([lm(1) lm(2)].*100)./100);
end;

%%
hits.powspctrm = (hits.powspctrm - M1);
misses.powspctrm = (misses.powspctrm - M2);

dum1 = squeeze(mean([hits.powspctrm],1));
dum2 = squeeze(mean([misses.powspctrm],1));

lm = [];
lm(1) = min(min( [dum1;dum2] ));
lm(2) = max(max( [dum1;dum2] ));

%%

figure;
for it = 1:24
    subplot(24/8,8,it);
    hold on;
    plot(misses.freq,squeeze(mean(mean(misses.powspctrm(:,it,:),2),1)),'k');
    plot(hits.freq,squeeze(mean(mean(hits.powspctrm(:,it,:),2),1)),'r');
    axis tight;
    ylim([lm]);
    set(gca,'YTick',[lm(1) lm(2)]);
    set(gca,'YTickLabel',round([lm(1) lm(2)].*100)./100);
end;

%%
hits.powspctrm = (hits.powspctrm)./M1;
misses.powspctrm = (misses.powspctrm)./M2;

dum1 = squeeze(mean([hits.powspctrm],1));
dum2 = squeeze(mean([misses.powspctrm],1));

lm = [];
lm(1) = min(min( [dum1;dum2] ));
lm(2) = max(max( [dum1;dum2] ));

%%
figure;
for it = 1:24
    subplot(24/8,8,it);
    hold on;
    plot(misses.freq,squeeze(mean(misses.powspctrm(:,it,:),1)),'k');
    plot(hits.freq,squeeze(mean(hits.powspctrm(:,it,:),1)),'r');
    axis tight;
    ylim([lm]);
    set(gca,'YTick',[lm(1) lm(2)]);
    set(gca,'YTickLabel',round([lm(1) lm(2)].*100)./100);
end;

%%
cfg                     = [];
cfg.parameter           = 'powspctrm';

[dum1] = ft_appendfreq( cfg , pow1{1}{1} , pow1{2}{1} , pow1{3}{1} );
[dum2] = ft_appendfreq( cfg , pow2{1}{1} , pow2{2}{1} , pow2{3}{1} );

%%
nchans = 16;%length(dum1.label);
figure;
for it = 1:16%length(dum1.label)
    
    y1 = squeeze(mean(dum1.powspctrm(:,it,:),1));
    y2 = squeeze(mean(dum2.powspctrm(:,it,:),1));
     
    subplot(nchans/8,8,it);
    hold on;
    plot(dum1.freq,(y1-y2)./(y1+y2));
    %plot(dum1.freq,y2,'r');
    axis tight;
    set(gca,'XTick',[1 10 20 30]);
end;