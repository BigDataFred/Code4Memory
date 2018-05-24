  
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