function [dat] = makeAVG4BF(dat)

%%
chanLabel = cell( 1 , length( dat.label ) );
for it = 1:length( dat.label )
    x = dat.label{it};
    chanLabel(it) = {x(1:regexp(x,'\d{1,2}')-1)};
end;

%%
BFid = unique(chanLabel);

%%
dum = cell( 1 , length(BFid) );
for it = 1:length( BFid )

    selIx = find(strcmp(chanLabel,BFid(it)));
    
    cfg                     = [];
    cfg.channel             = dat.label(selIx);
    cfg.avgoverchan         = 'yes';
    
    [dum{it}] = ft_selectdata( cfg , dat );
    
end;

[dat] = ft_appenddata( [], dum{:} );
dat.label = BFid;

return;