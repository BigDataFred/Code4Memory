%%
for it = 1:length( TFR )
    
    cfg                     = [];
    cfg.frequency           = [min(TFR{it}.freq) 19];
    cfg.latency             = [-1 max(TFR{it}.time)];
    
    [TFRlow] = ft_selectdata( cfg , TFR{it} );
    
    cfg                     = [];
    cfg.frequency           = [20 max(TFR{it}.freq)];
    cfg.latency             = [-1 max(TFR{it}.time)];
    
    [TFRhigh] = ft_selectdata( cfg , TFR{it} );
    
    %%
    
    %             %%
    %             cfg                     = [];
    %             cfg.baseline            = [-1 -0.05];
    %             cfg.baselinetype        = 'relchange';
    %
    %             [TFRlow] = ft_freqbaseline( cfg , TFRlow);
    %             [TFRhigh] = ft_freqbaseline( cfg , TFRhigh);
    
    %%
    cfg                     = [];
    cfg.frequency           = [55 65];
    cfg.avgoverfreq         = 'yes';
    
    [powH] = ft_selectdata( cfg , TFRhigh );
    
    %%
    figure;
    subplot(5,1,1:3);
    hold on;
    imagesc(TFRhigh.time,TFRhigh.freq,squeeze(mean((TFRhigh.powspctrm),1)));
    axis xy;
    plot([0 0],[min(TFRhigh.freq) max(TFRhigh.freq)],'w');
    plot([2 2],[min(TFRhigh.freq) max(TFRhigh.freq)],'w');
    axis tight;
    chName = TFRhigh.label;
    chName{:}(regexp(chName{:},'_')) = [];
    
    title([chName{:},' (',num2str(size(TFRhigh.powspctrm,1)),' trials)']);
    
    subplot(5,1,4:5);
    hold on;
    imagesc(TFRlow.time,TFRlow.freq,squeeze(mean((TFRlow.powspctrm),1)));
    axis xy;
    plot([0 0],[min(TFRlow.freq) max(TFRlow.freq)],'w');
    plot([2 2],[min(TFRlow.freq) max(TFRlow.freq)],'w');
    axis tight;
    
    %         %%
    %         M = squeeze(mean(powH.powspctrm,1));
    %         SE = std(powH.powspctrm,0,1)./sqrt(size(powH.powspctrm,1)-1);
    %
    %         figure;
    %         errorbar(powH.time,M,SE,'b');
    %         axis tight;
    %         title(powH.label);
    
    %%
    cfg                     = [];
    cfg.avgovertime         = 'yes';
    cfg.nanmean             = 'yes';
    
    [dum1]  = ft_selectdata( cfg , TFRlow );
    
    cfg                     = [];
    cfg.avgovertime         = 'yes';
    cfg.nanmean             = 'yes';
    
    [dum2]  = ft_selectdata( cfg , TFRhigh );
    
    cfg                     = [];
    cfg.variance            = 'yes';
    cfg.jackknife           = 'yes';
    
    [dum1]  = ft_freqdescriptives( cfg , dum1 );
    [dum2]  = ft_freqdescriptives( cfg , dum2 );
    
    %%
    figure;
    subplot(121);
    hold on;
    plot(dum1.freq,dum1.powspctrm);
    plot(dum1.freq,dum1.powspctrm-dum1.powspctrmsem);
    plot(dum1.freq,dum1.powspctrm+dum1.powspctrmsem);
    axis tight;
    title('Low frequencies');
    
    subplot(122);
    hold on;
    plot(dum2.freq,dum2.powspctrm);
    plot(dum2.freq,dum2.powspctrm-dum2.powspctrmsem);
    plot(dum2.freq,dum2.powspctrm+dum2.powspctrmsem);
    axis tight;
    title('High frequencies');
    
end;

%%
for it = length( TFR )
    
    %%
    x = seshTrl{it}';
    ix ={};
    for jt = 1:length( x )
        if jt ==1
            ix{jt} = 1:x(jt);
        else
            if ~isempty( ix{jt-1})
                ix{jt} = ix{jt-1}(end)+1:ix{jt-1}(end)+x(jt);
            else
                ix{jt} = ix{jt-2}(end)+1:ix{jt-2}(end)+x(jt);
            end;
        end;
    end;
    
    chck = zeros(1,length(ix));
    for jt = 1:length(ix)
        if isempty(ix(jt))
            chck(jt) = 1;
        end;
    end;

    ix(find(chck==1)) = [];
    
    %%
    for jt = 1:length(ix)
        
        cfg                     = [];
        cfg.frequency           = [min(TFR{it}.freq) 19];
        cfg.latency             = [-1 max(TFR{it}.time)];
        cfg.trials              = ix{jt};
        
        [TFRlow] = ft_selectdata( cfg , TFR{it} );
        
        cfg                     = [];
        cfg.frequency           = [20 max(TFR{it}.freq)];
        cfg.latency             = [-1 max(TFR{it}.time)];
        cfg.trials              = ix{jt};
        
        [TFRhigh] = ft_selectdata( cfg , TFR{it} );
        
        %%
        
%             %%
%             cfg                     = [];
%             cfg.baseline            = [-1 -0.05];
%             cfg.baselinetype        = 'relchange';
%         
%             [TFRlow] = ft_freqbaseline( cfg , TFRlow);
%             [TFRhigh] = ft_freqbaseline( cfg , TFRhigh);
        
        %%
        cfg                     = [];
        cfg.frequency           = [55 65];
        cfg.avgoverfreq         = 'yes';
        
        [powH] = ft_selectdata( cfg , TFRhigh );
        
        %%
        figure;
        subplot(5,1,1:3);
        hold on;
        imagesc(TFRhigh.time,TFRhigh.freq,squeeze(mean((TFRhigh.powspctrm),1)));
        axis xy;
        plot([0 0],[min(TFRhigh.freq) max(TFRhigh.freq)],'w');
        plot([2 2],[min(TFRhigh.freq) max(TFRhigh.freq)],'w');
        axis tight;
        chName = TFRhigh.label;
    chName{:}(regexp(chName{:},'_')) = [];
    
    title([chName{:},' (',num2str(size(TFRhigh.powspctrm,1)),' trials)']);
        subplot(5,1,4:5);
        hold on;
        imagesc(TFRlow.time,TFRlow.freq,squeeze(mean((TFRlow.powspctrm),1)));
        axis xy;
        plot([0 0],[min(TFRlow.freq) max(TFRlow.freq)],'w');
        plot([2 2],[min(TFRlow.freq) max(TFRlow.freq)],'w');
        axis tight;
        
%         %%
%         M = squeeze(mean(powH.powspctrm,1));
%         SE = std(powH.powspctrm,0,1)./sqrt(size(powH.powspctrm,1)-1);
%         
%         figure;
%         errorbar(powH.time,M,SE,'b');
%         axis tight;
%         title(powH.label);
        
        %%
        cfg                     = [];
        cfg.avgovertime         = 'yes';
        cfg.nanmean             = 'yes';
        
        [dum1]  = ft_selectdata( cfg , TFRlow );
        
        cfg                     = [];
        cfg.avgovertime         = 'yes';
        cfg.nanmean             = 'yes';
        
        [dum2]  = ft_selectdata( cfg , TFRhigh );
        
        cfg                     = [];
        cfg.variance            = 'yes';
        cfg.jackknife           = 'yes';
        
        [dum1]  = ft_freqdescriptives( cfg , dum1 );
        [dum2]  = ft_freqdescriptives( cfg , dum2 );
        
        %%
        figure;
        subplot(121);
        hold on;
        plot(dum1.freq,dum1.powspctrm);
        plot(dum1.freq,dum1.powspctrm-dum1.powspctrmsem);
        plot(dum1.freq,dum1.powspctrm+dum1.powspctrmsem);
        axis tight;
        title('Low frequencies');
        
        subplot(122);
        hold on;
        plot(dum2.freq,dum2.powspctrm);
        plot(dum2.freq,dum2.powspctrm-dum2.powspctrmsem);
        plot(dum2.freq,dum2.powspctrm+dum2.powspctrmsem);
        axis tight;
        title('High frequencies');
        
    end;
    
end;

%%
[raw] = eval(micLab{end});

cfg                     = [];
cfg.trials              = ix{end};
cfg.latency             = [-1 0];

[dum] = ft_selectdata( cfg, raw );

dt = 350;

cfg                     = [];
cfg.bpfilter            ='yes';
cfg.bpfreq              = [50 70];

[filt] = ft_preprocessing( cfg , dum );

epch1 = [];
c = 0;
for it = 1:length( filt.trial )
    fprintf([num2str(it),'/',num2str(length(filt.trial))]);
    
    [mix1] = local_max(filt.trial{it});
    [mix2] = local_max(-filt.trial{it});
    
    for jt = 1:length(mix1)
        selIx = mix1(it)-dt:mix1(it)+dt;
        
        if (min(selIx)>0) && (max(selIx) < length(dum.trial{it}))
            c = c + 1;
            epch1(c,:) = dum.trial{it}(selIx);
        end;
        
    end;
    
    fprintf('\n');
end;


cfg                     = [];
cfg.trials              = ix{end};
cfg.latency             = [2 max(raw.time{1})];
[dum] = ft_selectdata( cfg, raw );

dt = 350;

cfg                     = [];
cfg.bpfilter            ='yes';
cfg.bpfreq              = [50 70];

[filt] = ft_preprocessing( cfg , dum );

epch2 = [];
c = 0;
for it = 1:length( filt.trial )
    fprintf([num2str(it),'/',num2str(length(filt.trial))]);
    
    [mix1] = local_max(filt.trial{it});
    [mix2] = local_max(-filt.trial{it});
    
    for jt = 1:length(mix1)
        selIx = mix1(it)-dt:mix1(it)+dt;
        
        if (min(selIx)>0) && (max(selIx) < length(dum.trial{it}))
            c = c + 1;
            epch2(c,:) = dum.trial{it}(selIx);
        end;
        
    end;
    
    fprintf('\n');
end;

%%
dum1 = struct;
dum1.label = {'sta_chan_base'};
dum1.fsample = 1e3;
dum1.cfg = [];
for it = 1:size( epch1 ,1)
    dum1.time{it} = [-dt:dt];
    dum1.time{it} = dum1.time{it}.*1e-3;
    dum1.trial{it} = epch1(it,:);
end;

dum2 = struct;
dum2.label = {'sta_chan_post'};
dum2.fsample = 1e3;
dum2.cfg = [];
for it = 1:size( epch2 ,1)
    dum2.time{it} = [-dt:dt];
    dum2.time{it} = dum2.time{it}.*1e-3;
    dum2.trial{it} = epch2(it,:);
end;


%%
cfg                     = [];
cfg.method              = 'mtmfft';
cfg.pad                 = 'maxperlen';
cfg.output              = 'pow';
cfg.tapsmofrq           = 1.5;

[pow1] = ft_freqanalysis( cfg , dum1 );

cfg                     = [];
cfg.method              = 'mtmfft';
cfg.pad                 = 'maxperlen';
cfg.output              = 'pow';
cfg.tapsmofrq           = 1.5;

[pow2] = ft_freqanalysis( cfg , dum2 );