%%
restoredefaultpath;
addpath('~rouxf/fieldtrip-20161009/');
addpath(genpath('~rouxf/AnalysisFred/EM/'));
ft_defaults;
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
TFR = {};  
for it = 1:length(d)
    
    f = dir([p2d,d(it).name , filesep , 'lfp_dat', filesep , '*_EMtask_preprocLFPdata_ds.mat']);
    load([ p2d , d(it).name , filesep , 'lfp_dat', filesep , f.name ]);
    
    cfg = [];
    cfg.latency = [-2 3+min(y)];
    
    [lfp_data] = ft_selectdata( cfg , lfp_data );
    
    [TFR{it}] = time_freq_EM( lfp_data, 'low' , [11 13] );

end;

%%
cfg = [];
cfg.parameter = 'powspctrm';

[TFR] = ft_appendfreq( cfg, TFR{1}{1}, TFR{2}{1}, TFR{3}{1});
%TFR = TFR{1}{1};

cfg = [];
cfg.latency = [-.75 -.25];

[base] = ft_selectdata( cfg , TFR );

M = squeeze(mean(mean(base.powspctrm,4),1));
M = repmat(M,[1 1 length(TFR.time)]);

SD = squeeze(std(mean(base.powspctrm,4),1));
SD = repmat(SD,[1 1 length(TFR.time)]);

cfg = [];

[TFR] = ft_freqdescriptives( cfg , TFR );
%%
%Y1 = squeeze((TFR.powspctrm(1,:,:)-M(1,:,:))./SD(1,:,:));
%Y2 = squeeze((TFR.powspctrm(2,:,:)-M(2,:,:))./SD(2,:,:));
%Y1 = squeeze(mean((TFR.powspctrm-M)./SD,1));
Y1 = squeeze(mean((TFR.powspctrm),1));
%Y2 = squeeze((TFR.powspctrm(2,:,:)));

Y1 = (Y1-min(min(Y1)))./(max(max(Y1)) - min(min(Y1)));
%Y2 = (Y2-min(min(Y2)))./(max(max(Y2)) - min(min(Y2)));

%%
load('/home/rouxf/AnalysisFred/EM/colormapNOV16.mat');

figure;
%subplot(2,1,1);
a = gca;
pcolor(TFR.time, TFR.freq,Y1);shading interp;lighting phong;
axis xy;
% subplot(2,1,2);
% a = [a gca];
% pcolor(TFR.time, TFR.freq,Y2);shading interp;lighting phong;
% axis xy;

for it = 1:length(a)
    hold(a(it),'on');
    plot(a(it),[0 0],[min(TFR.freq) max(TFR.freq)],'w');
    plot(a(it),[2 2],[min(TFR.freq) max(TFR.freq)],'w');
    colormap(a(it),MAP);
end;
set(a,'Xlim',[-.75 4]);
set(gca,'XTick',[-.75 0 2 4]);
set(gcf,'Color','w');
cb = colorbar;
set(cb,'YTick',[0 .5 1]);
