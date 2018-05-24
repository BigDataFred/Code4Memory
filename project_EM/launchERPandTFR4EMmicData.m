%% set the path environment
addpath('/home/rouxf/tbx/fieldtrip-20170618/');
ft_defaults;

addpath('/home/rouxf/prj/Bham/code/mcode/project_EM/');
addpath('/home/rouxf/prj/Bham/code/mcode/utils/LFP/');
addpath('/home/rouxf/prj/Bham/code/mcode/visualize/continuous/');

%%
rpath = '/media/rouxf/rds-share/iEEG_DATA/MICRO/P08/fVSpEM/';

sesh = dir(rpath);
sesh(1:2) = [];

x = cell(1,length(sesh));
for it = 1:length(sesh)
    x(it) = {sesh(it).name};
end;
sesh = x;

%%
micDat = {};
micLab = {};

c = 0;
for it = 1:length( sesh )
    
    fprintf([sesh{it}]);
    
    spath = [rpath,sesh{it},filesep,'lfp_dat',filesep];
    lpath = [rpath,sesh{it},filesep,'log_dat',filesep];
    
    [micFiles] = dir([spath,'*_preprocANDrejectedIED_stimlocked.mat']);
    [logFiles] = dir([lpath,'*_LogFile_EMtask_LogDat.mat']);
    
    if ~isempty(logFiles)
        
        c =c+1;
        dat = load([lpath,logFiles(1).name]);
        [trl_all] = extract_Enc_trl(dat.LogDat1);
        
        for jt = 1:length( micFiles )
            micLab(jt) = {micFiles(jt).name(regexp(micFiles(jt).name,'CSC')+4:regexp(micFiles(jt).name,'201\d{1}')-2)};
        end;
        
        for jt = 1:length( micLab )
            
            dat = load([spath,micFiles(jt).name]);
            micDat{jt,c} = dat.LFPdatMW;
            selIx = dat.selIx;
            
            dum = intersect(trl_all,selIx);
            selIx = find(ismember(selIx,dum));
            
            if ~isempty(selIx)
                cfg                     = [];
                cfg.trials              = sort(selIx);
                
                [micDat{jt,c}] = ft_selectdata( cfg , micDat{jt,c} );
            else
                [micDat{jt,c}] = {};
            end;            
            
        end;
        
        clear dat;
    else
        ['toto',sesh{it}]
    end;
    
    fprintf('\n');
    
end;

%%
seshTrl = [];
for jt = 1:length( micLab )
    x = micDat(jt,:);
    
    chck = zeros(1,length(x));
    for it = 1:length(x)
        if isempty(x{it})
            chck(it) = 1;
        end;
    end;
    x(chck==1) = [];
    
    eval([micLab{jt},' = ft_appenddata( [], x{:} )']);
    
    dum= [];
    for it = 1:length(x)
        dum = [dum length(x{it}.trial)];
    end;
    seshTrl{jt} = dum;
    
end;

%%
erpDat = {};
TFR = {};
for jt = 1:length( micLab )
    
    [erpDat{jt},TFR{jt}] = compute_ERP_and_TFR_macroData_EM(eval(micLab{jt}));
    
    cfg                     = [];
    cfg.output              = 'pow';
    cfg.keeptrials          = 'yes';
    cfg.toi                 = [-1 4.475];
    
    [TFR{jt},slp]=sh_subtr1of(cfg,TFR{jt});
    
end;

%%

