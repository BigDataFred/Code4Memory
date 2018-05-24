function [MWdat,IEDdat,chID,origTime] = extractMWandIEDdat(micDat,micDatMode,MicTrlIx1,MicTrlIx2,IEDidx)

MWdat = cell(1,length(micDat));
IEDdat = cell(1,length(micDat));
for it = 1:length( micDat )
    
    if isfield(micDat{it},'LFPdatMW') && micDatMode == 2        
        cfg                     = [];
        cfg.trials              = MicTrlIx2;
        [MWdat{it}] = ft_selectdata( cfg , micDat{it}.LFPdatMW );          
    else
        cfg                     = [];
        cfg.trials              = MicTrlIx1;
        [MWdat{it}] = ft_selectdata( cfg , micDat{it}.save_data{1}{1}{1} );
        cfg                     = [];
        cfg.trials              = IEDidx;
        [IEDdat{it}] = ft_selectdata( cfg , micDat{it}.save_data{1}{1}{1} );
    end;
    
end;

if length(MWdat)>1
    MWdat = ft_appenddata([], MWdat{:});
    try
        IEDdat = ft_appenddata([], IEDdat{:});
    catch
        IEDdat = [];
    end;
else
    MWdat  = MWdat{1};
    IEDdat = IEDdat{1};
end;


chID = {};
for it = 1:length(MWdat.label)
    ix = 1:regexp(MWdat.label{it},'\d{1,2}')-1;
    chID(it) = {MWdat.label{it}(ix)};
end;

origTime = MWdat.time{1};
