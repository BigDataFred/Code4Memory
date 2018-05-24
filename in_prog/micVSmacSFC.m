function micVSmacSFC

%%
addpath('/home/rouxf/tbx/fieldtrip-20170618/');
ft_defaults;

addpath('/home/rouxf/prj/Bham/code/mcode/project_EM/');
addpath('/home/rouxf/prj/Bham/code/mcode/utils/LFP/');
addpath('/home/rouxf/prj/Bham/code/mcode/visualize/continuous/');

%%
[p2MIC] = '/media/rouxf/rds-share/iEEG_DATA/MICRO/P08/fVSpEM/2017-05-31_10-52-16/lfp_dat/';
[micFiles] = dir([p2MIC,'*stimlocked.mat']);

%[p2MAC] = '/media/rouxf/rds-share/iEEG_DATA/MACRO/P02/fVSpEM/2016-07-17_13-45-06/';
%[macFile] = dir([p2MAC,'*_PreprocMacroLFP_fVSpEM_*.mat']);

mode = 0;% 0 = micro-data only (default), 1 = macro-data only, 2 = both micro & macro

%% load & preprocess the MW-data
if ( mode == 0 ) || ( mode == 2 )
    
    [micDat] = cell( 1 , length(micFiles) );
    for it = 1:length(micFiles)
        
        fprintf([num2str(it),'/',num2str( length(micFiles) )]);
        [dat] = load([p2MIC,micFiles(it).name]);
        [micDat{it}] = dat.save_data{1}{1}{1};
        fprintf('\n');
        
    end;
    
    [micDat] = ft_appenddata( [] , micDat{:} );
    
    %% average the signals of the MWs of each BF-shank
    [micDat] = makeAVG4BF(micDat);
    
    %% preproc Micro-data
    cfg                     = [];
    cfg.detrend             = 'yes';
    cfg.demean              = 'yes';
    cfg.bpfilter            = 'yes';
    cfg.bpfilttype          = 'fir';
    cfg.bpfreq              = [0.5 200];
    cfg.padtype             = 'mirror';
    cfg.padding             = 5;
    
    [micDat] = ft_preprocessing( cfg , micDat );
    
end;

%% load & preprocess the Macro-data
if ( mode == 1 ) || ( mode == 2 )
    
    [dat] = load([p2MAC,macFile.name]);
    [macDat] = dat.macroLFPdat;
    
    %%
    cfg                     = [];
    cfg.channel             = {'*H*'};
    
    [macDat] = ft_selectdata( cfg , macDat );
    
    %%
    [macDat] = extractDat4DistalMacContact(macDat);
    
    %%
    cfg                     = [];
    cfg.detrend             = 'yes';
    cfg.demean              = 'yes';
    cfg.bpfilter            = 'yes';
    cfg.bpfilttype          = 'fir';
    cfg.bpfreq              = [0.5 200];
    cfg.padtype             = 'mirror';
    cfg.padding             = 5;
    
    [macDat] = ft_preprocessing( cfg , macDat );    
    
end;

%%
nsamp = length(micDat.trial{1});
ntrl = length( micDat.trial );

samp = zeros(ntrl,2);
for it = 1:ntrl
    samp(it,1) = (it-1)*nsamp+1;
    samp(it,2) = (it-1)*nsamp+nsamp;
end;

%%
[~,trlIx] = removeIEDtrl(micDat,.25,150,'n');

%%
chanDat = cell( 1 , length(micDat.label) );
delIx = {};
for it = 1:length( micDat.label)
    
    delIx{it} = trlIx{it};
    
    selIx = setdiff(1:length(micDat.trial),trlIx{it});
    
    cfg                     = [];
    cfg.channel             = micDat.label;
    cfg.trials              = selIx;
    
    [chanDat{it}] = ft_selectdata( cfg , micDat );
%     
%     f = 0;
%     while f <1
%         
%         cfg                     = [];
%         cfg.method              = 'summary';
%         
%         [dum] = ft_rejectvisual( cfg , chanDat{it} );
%         
%         cfg                     = [];
%         cfg.method              = 'trial';
%         
%         [dum] = ft_rejectvisual( cfg , dum );
%         
%         [s] = input('Do you wish to keep the current selection [y/n]','s');
%         
%         if strcmp(s,'y')
%             f = 1;
%             chanDat{it} = dum;
%         end;
%         
%     end;
%     
%     if isfield(dum.cfg.artfctdef,'summary')
%         x = dum.cfg.artfctdef.summary.artifact;
%         ix = zeros(1,size(x,1));
%         
%         for jt = 1:size(x,1)
%             [dum,~] = find(ismember(samp,x(jt,:)));
%             ix(jt) = unique(dum);
%         end;
%         
%         delIx{it} = [delIx{it},ix];
%     end;

end;

%%
savepath = '';
savename = '';

save(savepath,savename,'micDat','delIx');

%%
return;

%%
function [chanIx,trlIx] = removeIEDtrl(dat,pctTrsh,dt,plt)

ntrl = length(dat.trial);

[dat,ixON,~] = runIEDdetection(dat,dt,plt);

[pct,trlChck] = calculateIEDPct(ntrl,ixON,pctTrsh,plt);

delIx = find( pct >= pctTrsh );
[chanIx] = setdiff(1:length(dat.label),delIx);
if isempty(chanIx);chanIx = delIx;end;

[trlIx] = trlChck;

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
