%%
ft = dir('/home/rouxf/tbx/fieldtrip-********');
addpath(['/home/rouxf/tbx/',ft.name]);
ft_defaults;

%%
p2d = '/media/rouxf/rds-share/iEEG_DATA/MACRO/P04/fVSpEM/';
load([p2d,'P04_macroData_fVSpEM_TFRandERP_2016-10-18_15-30-26.mat']);

%%
TFR = save_data{1}{1};
ERP = save_data{1}{2};

%%
cfg = [];
cfg.channel = {'LEn*' 'LmH*' 'RmH*'  'Lwp-i*'};

[TFR] = ft_selectdata( cfg , TFR );
[ERP] = ft_selectdata( cfg , ERP );

%%
% cfg = [];
% cfb.baselinetype ='absolute';
% cfg.baseline = [-1 0];
% 
% [TFR] = ft_freqbaseline(cfg,TFR);

%%
cfg = [];
cfg.frequency = [min(TFR.freq) 19];

[lfPOW] = ft_selectdata(cfg,TFR);

cfg = [];
cfg.frequency = [20 max(TFR.freq)];

[hfPOW] = ft_selectdata(cfg,TFR);

%%
itc = zeros(size(ERP.trial));
for it = 1:length( ERP.label )
    
    dum = squeeze( ERP.trial(:,it,:) );
    dum = hilbert(dum')';
    
    itc(:,it,:) = dum;
    
end;

itc = itc./abs(itc);
itc = abs(sum(itc,1))./size(itc,1);
itc = squeeze( itc );

%%
figure;
for it = 1:size(ERP.trial,2);
    it
    subplot(12,1,1:3);
    hold on;
    imagesc(ERP.time,1:size(ERP.trial,1),squeeze(ERP.trial(:,it,:)));
    %ca = caxis;
    plot([0 0],[1 size(ERP.trial,1)],'w');
    plot([2 2],[1 size(ERP.trial,1)],'w');
    axis tight;
    
    title(ERP.label(it));
    
    subplot(12,1,4);
    [ax,h1,h2] = plotyy(ERP.time,ERP.avg(it,:),ERP.time,itc(it,:));
    axis(ax,'tight');
    
    subplot(12,1,6:8);
    hold on;
    imagesc(hfPOW.time,hfPOW.freq,squeeze(mean(mean((hfPOW.powspctrm(:,it,:,:)),2),1)));
    plot([0 0],[min(hfPOW.freq) max(hfPOW.freq)],'w');
    plot([2 2],[min(hfPOW.freq) max(hfPOW.freq)],'w');
    axis xy;axis tight;ylim([min(hfPOW.freq) max(max(hfPOW.freq))]);
    %caxis([-50 50]);
    
    subplot(12,1,10:12);
    hold on;
    imagesc(lfPOW.time,lfPOW.freq,squeeze(mean(mean((lfPOW.powspctrm(:,it,:,:)),2),1)));
    plot([0 0],[min(lfPOW.freq) max(lfPOW.freq)],'w');
    plot([2 2],[min(lfPOW.freq) max(lfPOW.freq)],'w');
    axis xy;axis tight;ylim([min(lfPOW.freq) max(max(lfPOW.freq))]);
    %caxis([-250 250]);
    
    pause;
    clf;
end;

%%
%[erpDat,TFR] = compute_ERP_and_TFR_macroData_EM(macDat);

% %%
% itc = zeros(length(erpDat.label),length(erpDat.time));
% for it = 1:length(erpDat.label)
%     
%     x = squeeze( erpDat.trial(:,it,:) );
%     x = hilbert(x')';
%     x = x./abs(x);
%     itc(it,:) = abs(sum(x,1))./size(x,1);
% end;


% %%
% figure;
% for it = 1:length( erpDat.label)
%     
%     subplot(4,1,1:3);
%     hold on;
%     imagesc(erpDat.time,1:size(erpDat.trial,1),squeeze( erpDat.trial(:,it,:) ) );
%     plot([0 0],[1 size(erpDat.trial,1)],'w');
%     plot([2 2],[1 size(erpDat.trial,1)],'w');
%     axis xy;axis tight;
%     caxis([-100 100]);
%      
%     subplot(4,1,4);
%     hold on;
%     [ax,h1,h2] = plotyy(erpDat.time,erpDat.avg(it,:),erpDat.time,itc(it,:));
%     plot(ax(1),[0 0],[min(erpDat.avg(it,:)) max(erpDat.avg(it,:))],'k');
%     plot(ax(1),[2 2],[min(erpDat.avg(it,:)) max(erpDat.avg(it,:))],'k');
%     axis(ax,'tight');
%     
%     pause;
%     clf;
% end;
% 
% %%
% selIx1 = find(TFR.freq<19);
% selIx2 = find(TFR.freq>=19);
% 
% for it = 1:length(TFR.label)
%     
%     figure;
%     for jt = 1:size(TFR.powspctrm,1)
%         subplot(4,1,1:2);
%         imagesc(TFR.time,TFR.freq(selIx2),squeeze(mean(20*log10(TFR.powspctrm(jt,it,selIx2,:)),2)));
%         axis xy;
%         
%         subplot(4,1,3:4);
%         imagesc(TFR.time,TFR.freq(selIx1),squeeze(mean(20*log10(TFR.powspctrm(jt,it,selIx1,:)),2)));
%         axis xy;
%         
%         pause;
%         clf;
%         
%     end;
% end;
