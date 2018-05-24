%%
restoredefaultpath;
addpath('~rouxf/fieldtrip-20161009/');
ft_defaults;
%%
pID = 'P04';
p2f = [filesep,'media',filesep,'rouxf',filesep,'rds-share',filesep,'iEEG_DATA',filesep,'MICRO',filesep,pID,filesep,'fvsPEM',filesep];
d = dir([p2f,'2016*']);
%%
LFP_data = cell(1,length(d));
for it = 1:length(d)
    
    p2f2 = [p2f,d(it).name,filesep,'preproc',filesep];
    
    fn = dir([p2f2,'*_EMtask_preprocLFPdata.mat']);
    fn = fn.name;
    
    dat = load([p2f2,fn]);
    
    cfg                     = [];
    cfg.lpfilter            = 'yes';
    cfg.lpfreq              = 30;
    
    [dum] = ft_preprocessing(cfg , dat.LFP_data );
    clear dat;
    
    cfg                     = [];
    cfg.channel             = {'L*'}; 
    
    [dum] = ft_selectdata(cfg , dum );
    
    cfg                     = [];
    cfg.resamplefs          = 1200;
    
    [LFP_data{it}] = ft_resampledata(cfg , dum );
        
end;
%%
ERP = cell(1,length(d));
for it = 1:length(LFP_data)
    
    cfg                         = [];
    cfg.keeptrials              = 'yes';
    
    dum = ft_timelockanalysis(cfg , LFP_data{it} );
    
    RMS = sqrt(mean(abs(dum.trial).^2,3));
    
    trl_ix = [];
    for jt = 1:size(RMS,2)
        z = zscore(RMS(:,jt));
        trl_ix = [trl_ix;find(abs(z) >=2.5)];
    end;
    trl_ix = unique(trl_ix);
    
    cfg                     = []; 
    cfg.trials              = setdiff(1:size(dum.trial,1),trl_ix);
    
    [dum] = ft_selectdata(cfg , LFP_data{it} );
    
    cfg                     = [];
    cfg.keeptrials          = 'yes';
    
    [ERP{it}] = ft_timelockanalysis(cfg , dum );

end;

%%
sel_ix = cell(1,length(ERP));
pow = cell(1,length(ERP));

for it = 1:length(ERP)
    Phi = zeros(size(ERP{it}.trial));
    for jt = 1:size(ERP{it}.trial,1)
        for kt = 1:length(ERP{it}.label)
            Phi(jt,kt,:) =  hilbert(squeeze(ERP{it}.trial(jt,kt,:)));
        end;
    end;
    
    Phi = Phi(:,:,find(ERP{it}.time>=0));
    
    ITC = sum(Phi,1)./sum(abs(Phi),1);
    ITC = squeeze(abs(ITC));
    
    PLV = squeeze(1/size(Phi,1)*abs(sum(exp(1i.*angle(Phi)),1)));
    
    [v,ix] = max(ITC,[],2);
    ix = ix(find(v >=.5));
    v = v(find(v >=.5));
    
    [v,s_ix] = sort(v);
    
    ix = ix(s_ix);
    
    [ix1,ix2] = ind2sub(size(ITC),ix);
    sel_ix{it} = sort(ix1);
    
    cfg                 = [];
    cfg.latency         = [0 max(ERP{it}.time)];
    
    [dum] = ft_selectdata(cfg , ERP{it} );
        
    cfg                 = [];
    cfg.method          = 'mtmfft';
    cfg.pad             ='maxperlen';
    cfg.taper           ='dpss';
    cfg.tapsmofrq       =1;
    
    [pow{it}] = ft_freqanalysis(cfg, dum );
end;

%%

for it = 1:length(ERP)
    %%
    figure;    
    for jt = 1:length(ERP{it}.label)
                
        subplot(length(ERP{it}.label)/8,8,jt);
        hold on;
        h = area([0 2],[max(max(ERP{it}.avg)) max(max(ERP{it}.avg))],min(min(ERP{it}.avg)));
        set(h,'FaceColor',[.75 .75 .75],'EdgeColor',[.75 .75 .57],'ShowBaseLine','off');
        plot(ERP{it}.time,ERP{it}.avg(jt,:));
        axis tight;
        xlim([-1 4]);
        title(ERP{it}.label(jt));
        ylim([min(min(ERP{it}.avg))-5 max(max(ERP{it}.avg))+5]);
        
        if ismember(jt,sel_ix{it})
            set(gca,'LineWidth',3);
            box on;
        end;
        
    end;
    
    %%
    figure;
    k1 = 1;k2 = 2;
    for jt = 1:length(ERP{it}.label)
        
        n = length(ERP{it}.label);%size(ERP{it}.trial,1);
        
        subplot(n/2,2,jt);%k1
        hold on;
        imagesc(ERP{it}.time,1:size(ERP{it}.trial,1),squeeze(ERP{it}.trial(:,jt,:)));%sel_ix{it}(jt)
        wf = squeeze(ERP{it}.avg(jt,:));%sel_ix{it}(jt)
        wf = wf./max(wf);
        wf = wf.*10;
        plot(ERP{it}.time,wf+floor(n/2),'w');
        axis xy;axis tight;
        xlim([-1 max(ERP{it}.time)]);
        k1 = k1+2;
        %title(ERP{it}.label(jt));%sel_ix{it}(jt)
        
%         subplot(n,2,k2);
%         plot(pow{it}.freq,pow{it}.powspctrm(jt,:));%sel_ix{it}(jt)
%         axis tight;xlim([0 30]);
%         k2 = k2+2;
%         %title(pow{it}.label(jt));%sel_ix{it}(jt))
    end;
    
end;


