function [delIx,selIx] = detectIED2trlIx(lfpDat,tIx)

%% extract BF labels and ids
ix = regexp(lfpDat.chanLab,'\d{1}');
ix = [ix{:}]-1;
[BFlab] = cell(length(ix),1);
[hemiLab] = cell(length(ix),1);

for it = 1:length(ix)
    BFlab(it) = { lfpDat.chanLab{it}(1:ix(it)) };
    hemiLab(it) = {BFlab{it}(end)};
end;
[BFid] = unique(BFlab);
[hemiID] = unique(hemiLab);

[nsmp,ntrl] = size(lfpDat.LFPseg{1}(tIx,:));

delIx = cell(length( hemiID ),1);
selIx = cell(length( hemiID ),1);

for it = 1:length( hemiID )
    
    selIx{it} = find(strcmp(hemiLab,hemiID(it)));
    
    [dum] = struct;
    dum.label = lfpDat.chanLab(selIx{it});
    dum.trial = cell(1,ntrl);
    dum.time = cell(1,ntrl);
    for kt = 1:ntrl
        dum.trial{kt} = zeros(length(selIx{it}),nsmp);
        for jt = 1:length( selIx{it})
            dum.trial{kt}(jt,:) = lfpDat.LFPseg{selIx{it}(jt)}(tIx,kt)';
            dum.time{kt} = lfpDat.trlTime(tIx);
        end;
    end;
    
    cfg                     = [];
    cfg.demean              = 'yes';
    cfg.detrend             = 'yes';
    cfg.bpfilter            = 'yes';
    cfg.bpfreq              = [0.5 120];
    
    [dum] = ft_preprocessing( cfg, dum );
    
    
    cfg                     = [];
    cfg.viewmode            = 'vertical';
    
    [cfg] = ft_databrowser(cfg,dum);
    
%     cfg = [];
%     cfg.continuous = 'no';
%     cfg.artfctdef.zvalue.interactive = 'yes';
%     cfg.artfctdef.zvalue.cutoff = 4;
%     cfg.artfctdef.zvalue.channel = dum.label;
%     
%     [cfg, artifact] = ft_artifact_zvalue(cfg, dum);
    
    x1 = cfg.artfctdef.visual.artifact;%artifact;%cfg.artfctdef.zvalue.artifact;%
    for jt = 1:ntrl
        x2 = (jt-1)*nsmp+1:nsmp*jt;  
        for kt = 1:size(x1,1)
            if (any(ismember(x1(kt,1):x1(kt,2),x2)))
                delIx{it} = [delIx{it} jt];
            end;
        end;
    end;
    delIx{it} = unique(delIx{it});
        
end;

%%
return;