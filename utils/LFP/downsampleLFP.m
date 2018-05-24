function [LFP,dsFs] = downsampleLFP(LFP,trlTime,dsFs)

for it = 1:length( LFP )
    
    dum = [];
    dum.label = {'dumChan'};
    dum.trial = cell(1,size(LFP{it},2));
    dum.time = cell(1,size(LFP{it},2));
    for jt = 1:size( LFP{it},2)
        dum.trial{jt} = LFP{it}(:,jt)';
        dum.time{jt} = trlTime;
    end;
    
    cfg                     = [];
    cfg.resamplefs          = dsFs;
    
    [dum] = ft_resampledata( cfg, dum );
    
    cfg                     = [];
    cfg.keeptrials          = 'yes';
    
    [dum] = ft_timelockanalysis( cfg , dum );
    
    LFP{it} = [];
    LFP{it} = squeeze(dum.trial)';
    
end;
