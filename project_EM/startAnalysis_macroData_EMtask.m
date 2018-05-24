function startAnalysis_macroData_EMtask(pID,expMode)
%%
restoredefaultpath;
addpath(genpath('/home/rouxf/prj/Bham/code/mcode/'));
ft = dir(['/home/rouxf/tbx/fieldtrip-********']);
addpath(['/home/rouxf/tbx/',ft.name]);ft_defaults;

%%
if nargin ==0
    pID = 'P02';
    expMode = 'fVSpEM';
end;

%%
%preproc_macro_data_EM(pID,expMode)

%%
p2d = ['/media/rouxf/rds-share/iEEG_DATA/MACRO/',pID,'/',expMode,'/'];
savepath = p2d;

c = 0;
sel = [];
sesh = dir(p2d);
for it = 1:length(sesh)
    if regexp(sesh(it).name,'\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}')
        c =c+1;
        sel(c) = it;
    end;
end;
sesh = sesh(sel);

%%
sel = [];
c = 0;
for it = 1:length( sesh )
    chck = dir([p2d,sesh(it).name]);    
    if length(chck) >2
        c = c+1;
        sel(c) = it;
    end;
end;
sesh = sesh(sel);

%%
if isempty(gcp('nocreate'))
    parpool(length(sesh),'SpmdEnabled',false);
end;

%%
ERP   = cell(1,length(sesh));
TFR   = cell(1,length(sesh));
NLXts = cell(1,length(sesh));

parfor it = 1:length(sesh)
%for it = 1:length(sesh)    
    dat = load([p2d,sesh(it).name,filesep,pID,'_PreprocMacroLFP_',expMode,'_',sesh(it).name,'.mat'])
    
    %% get the labels of the macro-electrodes
    [chanID]    = extractMacroChanID(dat.macroLFPdat.label);
    chanID_     = unique(chanID);
    nChans      = length(chanID_);
    
    [biRefLFPdat] = referenceDat(dat.macroLFPdat, chanID_, chanID , 'bipolar' );
    
    %% compute ERPs and TFRs
    [ERP,TFR] = compute_ERP_and_TFR_macroData_EM(biRefLFPdat);
    [NLXts]   = sesh(it).name;
    
    %%
    readme = {'TFR','ERP',NLXts};
    savename = [pID,'_macroData_',expMode,'_TFRandERP_',sesh(it).name,'.mat'];
    par_save([savepath,savename],[],{TFR,ERP,readme});
    
end;

%%
delete(gcp);

