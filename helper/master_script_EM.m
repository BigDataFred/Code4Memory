%% set the path def
restoredefaultpath;
addpath(genpath('/home/rouxf/prj/Bham/code/mcode/'));% custom code
addpath(genpath('/home/rouxf/tbx/releaseDec2015/'));%needed to read Nlx data
addpath(genpath('/home/rouxf/tbx/osort-v3-rel/'));%needed for spike detection & sorting
addpath(genpath('/home/rouxf/tbx/chronux_2_11/'));% needed for spectral analysis
addpath(genpath('/home/rouxf/tbx/wave_clus/'));%needed for spike detection & sorting
fn = dir('/home/rouxf/tbx/fieldtrip-*');
addpath(['/home/rouxf/tbx/',fn.name]);% spectral analysis, signal processing, spike detection
ft_defaults;

%% open pool of workers
if isempty(gcp('nocreate'))
    parpool(36,'SpmdEnabled',false);
end;

% %% set params
% params                  = [];
% params.pID              = 'P02';
% params.p2Nlxdata = ['/media/rouxf/rds-share/Archive/MICRO/',params.pID,'/fvSpEM/'];%path 2 Nlx data
% params.savepath = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',params.pID,'/fvSpEM/'];%path for saving
% params.mode         = 'stimlocked';%
% params.pre          = 5;% baseline time range
% params.post         = 10;% params.post stim time range
% params.makeLFP      = 1;
% params.makeSpikes   = 1;
% params.makeMUA      = 0;
% params.rpath = '/home/rouxf/';
% params.Nlxdat = dir(params.p2Nlxdata);
% params.Nlxdat(1:2) = [];
% params.Nlxdat = {params.Nlxdat.name}';
%     
% %%
% seshs = dir([params.p2Nlxdata,'2016*']);
% {seshs(:).name}'
% 
% %% launch parallel preproc
% for it =1:length(seshs)
%         
%     %%
%     params.sesh         = it;
%     nNCsf               = dir([params.p2Nlxdata,filesep,seshs(it).name,'/*.ncs']);
%     params.CSC_idx      = [ 1:length( nNCsf ) ];
%     
%     %%
%     try
%         launch_par_preproc_EM(params);
%     catch ME        
%         fid = fopen([params.savepath,'errorLog_',seshs(it).name,'.txt'],'w');
%         fprintf(fid,ME.message);
%         fclose(fid);
%     end;
% 
% end;
% 
% %% set params
% params                  = [];
% params.pID              = 'P04';
% params.p2Nlxdata = ['/media/rouxf/rds-share/Archive/MICRO/',params.pID,'/fvSpEM/'];%path 2 Nlx data
% params.savepath = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',params.pID,'/fvSpEM/'];%path for saving
% params.mode         = 'stimlocked';%
% params.pre          = 5;% baseline time range
% params.post         = 10;% params.post stim time range
% params.makeLFP      = 1;
% params.makeSpikes   = 1;
% params.makeMUA      = 0;
% params.rpath = '/home/rouxf/';
% params.Nlxdat = dir(params.p2Nlxdata);
% params.Nlxdat(1:2) = [];
% params.Nlxdat = {params.Nlxdat.name}';
% 
% %%
% seshs = dir([params.p2Nlxdata,'2016*']);
% {seshs(:).name}
% 
% %% launch parallel preproc
% for it = 1:length(seshs)
%     
%     %%
%     params.sesh         = it;
%     nNCsf               = dir([params.p2Nlxdata,filesep,seshs(it).name,'/*.ncs']);
%     params.CSC_idx      = [ 1:length( nNCsf ) ];
%     
%     %%
%     try
%         launch_par_preproc_EM(params);
%     catch ME
%         fid = fopen([params.savepath,'errorLog_',seshs(it).name,'.txt'],'w');
%         fprintf(fid,ME.message);
%         fclose(fid);
%     end;
% 
% end;
% 
% %% set params
% params                  = [];
% params.pID              = 'P05';
% params.p2Nlxdata = ['/media/rouxf/rds-share/Archive/MICRO/',params.pID,'/fvSpEM/'];%path 2 Nlx data
% params.savepath = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',params.pID,'/fvSpEM/'];%path for saving
% params.mode         = 'stimlocked';%
% params.pre          = 5;% baseline time range
% params.post         = 10;% params.post stim time range
% params.makeLFP      = 1;
% params.makeSpikes   = 1;
% params.makeMUA      = 0;
% params.rpath = '/home/rouxf/';
% params.Nlxdat = dir(params.p2Nlxdata);
% params.Nlxdat(1:2) = [];
% params.Nlxdat = {params.Nlxdat.name}';
% 
% %%
% seshs = dir([params.p2Nlxdata,'2016*']);
% {seshs(:).name}'
% 
% %% launch parallel preproc
% for it = 1:length(seshs)
%           
%     %%
%     params.sesh         = it;
%     nNCsf               = dir([params.p2Nlxdata,filesep,seshs(it).name,'/*.ncs']);
%     params.CSC_idx      = [ 1:length( nNCsf ) ];
%     
%     %%
%     try
%         launch_par_preproc_EM(params);
%     catch ME
%         fid = fopen([params.savepath,'errorLog_',seshs(it).name,'.txt'],'w');
%         fprintf(fid,ME.message);
%         fclose(fid);
%     end;
% 
% end;
%  
% %% set params
% params                  = [];
% params.pID              = 'P22AMS';
% params.p2Nlxdata = ['/media/rouxf/rds-share/Archive/MICRO/',params.pID,'/cnEM/'];%path 2 Nlx data
% params.savepath = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',params.pID,'/cnEM/'];%path for saving
% params.mode         = 'stimlocked';%'resplocked';%
% params.pre          = 5;% baseline time range
% params.post         = 10;% params.post stim time range
% params.makeLFP      = 1;
% params.makeSpikes   = 1;
% params.makeMUA      = 0;
% params.rpath = '/home/rouxf/';
% params.Nlxdat = dir(params.p2Nlxdata);
% params.Nlxdat(1:2) = [];
% params.Nlxdat = {params.Nlxdat.name}';
% 
% %%
% seshs = dir([params.p2Nlxdata,'2016*']);
% {seshs(:).name}'
% 
% %% launch parallel preproc
% for it = 1:length(seshs)
%     
%     %%
%     params.sesh         = it;
%     nNCsf               = dir([params.p2Nlxdata,filesep,seshs(it).name,'/*.ncs']);
%     params.CSC_idx      = [ 1:length( nNCsf ) ];
%     
%     %%
%     try
%         launch_par_preproc_EM(params);
%     catch ME
%         fid = fopen([params.savepath,'errorLog_',seshs(it).name,'.txt'],'w');
%         fprintf(fid,ME.message);
%         fclose(fid);
%     end;
% 
% end;
% 
% %% set params
% params                  = [];
% params.pID              = 'P23AMS';
% params.p2Nlxdata = ['/media/rouxf/rds-share/Archive/MICRO/',params.pID,'/cnEM/'];%path 2 Nlx data
% params.savepath = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',params.pID,'/cnEM/'];%path for saving
% params.mode         = 'stimlocked';%'resplocked';%
% params.pre          = 5;% baseline time range
% params.post         = 10;% params.post stim time range
% params.makeLFP      = 1;
% params.makeSpikes   = 1;
% params.makeMUA      = 0;
% params.rpath = '/home/rouxf/';
% params.Nlxdat = dir(params.p2Nlxdata);
% params.Nlxdat(1:2) = [];
% params.Nlxdat = {params.Nlxdat.name}';
% 
% %%
% seshs = dir([params.p2Nlxdata,'2016*']);
% {seshs(:).name}'
% 
% %% launch parallel preproc
% for it = 1:length(seshs)
%         
%     %%
%     params.sesh         = it;
%     nNCsf               = dir([params.p2Nlxdata,filesep,seshs(it).name,'/*.ncs']);
%     params.CSC_idx      = [ 1:length( nNCsf ) ];
%     
%     %%
%     try
%         launch_par_preproc_EM(params);
%     catch ME
%         fid = fopen([params.savepath,'errorLog_',seshs(it).name,'.txt'],'w');
%         fprintf(fid,ME.message);
%         fclose(fid);
%     end;
% 
% end;
% 
% %% set params
% params                  = [];
% params.pID              = 'P07';
% params.p2Nlxdata = ['/media/rouxf/rds-share/Archive/MICRO/',params.pID,'/fVSpEM/'];%path 2 Nlx data
% params.savepath = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',params.pID,'/fVSpEM/'];%path for saving
% params.mode         = 'stimlocked';%
% params.pre          = 5;% baseline time range
% params.post         = 10;% params.post stim time range
% params.makeLFP      = 1;
% params.makeSpikes   = 1;
% params.makeMUA      = 0;
% params.rpath = '/home/rouxf/';
% params.Nlxdat = dir(params.p2Nlxdata);
% params.Nlxdat(1:2) = [];
% params.Nlxdat = {params.Nlxdat.name}';
% 
% %%
% seshs = dir([params.p2Nlxdata,'2017*']);
% {seshs(:).name}'
% 
% %% launch parallel preproc
% for it = 1:length(seshs)
%         
%     %%
%     params.sesh         = it;
%     nNCsf               = dir([params.p2Nlxdata,filesep,seshs(it).name,'/*.ncs']);
%     params.CSC_idx      = [ 1:length( nNCsf ) ];
%     
%     %%
%     try
%         launch_par_preproc_EM(params);
%     catch ME
%         fid = fopen([params.savepath,'errorLog_',seshs(it).name,'.txt'],'w');
%         fprintf(fid,ME.message);
%         fclose(fid);
%     end;
% 
% end;
% 
% %% set params
% params                  = [];
% params.pID              = 'P07';
% params.p2Nlxdata = ['/media/rouxf/rds-share/Archive/MICRO/',params.pID,'/cnEM/'];%path 2 Nlx data
% params.savepath = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',params.pID,'/cnEM/'];%path for saving
% params.mode         = 'stimlocked';%
% params.pre          = 5;% baseline time range
% params.post         = 10;% params.post stim time range
% params.makeLFP      = 1;
% params.makeSpikes   = 1;
% params.makeMUA      = 0;
% params.rpath = '/home/rouxf/';
% params.Nlxdat = dir(params.p2Nlxdata);
% params.Nlxdat(1:2) = [];
% params.Nlxdat = {params.Nlxdat.name}';
% 
% %%
% seshs = dir([params.p2Nlxdata,'2017*']);
% {seshs(:).name}'
% 
% %% launch parallel preproc
% for it = 1:length(seshs)
%         
%     %%
%     params.sesh         = it;
%     nNCsf               = dir([params.p2Nlxdata,filesep,seshs(it).name,'/*.ncs']);
%     params.CSC_idx      = [ 1:length( nNCsf ) ];
%     
%     %%
%     try
%         launch_par_preproc_EM(params);
%     catch ME
%         fid = fopen([params.savepath,'errorLog_',seshs(it).name,'.txt'],'w');
%         fprintf(fid,ME.message);
%         fclose(fid);
%     end;
% 
% end;
% 
% %% set params
% params                  = [];
% params.pID              = 'P08';
% params.p2Nlxdata = ['/media/rouxf/rds-share/Archive/MICRO/',params.pID,'/cnEM/'];%path 2 Nlx data
% params.savepath =  ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',params.pID,'/cnEM/'];%path for saving
% params.mode         = 'stimlocked';%
% params.pre          = 5;% baseline time range
% params.post         = 10;% params.post stim time range
% params.makeLFP      = 1;
% params.makeSpikes   = 1;
% params.makeMUA      = 0;
% params.rpath = '/home/rouxf/';
% params.Nlxdat = dir(params.p2Nlxdata);
% params.Nlxdat(1:2) = [];
% params.Nlxdat = {params.Nlxdat.name}';
% 
% %%
% seshs = dir([params.p2Nlxdata,'2017*']);
% {seshs(:).name}'
% 
% %% launch parallel preproc
% for it = 1:length(seshs)
%             
%     %%
%     params.sesh         = it;
%     nNCsf               = dir([params.p2Nlxdata,filesep,seshs(it).name,'/*.ncs']);
%     params.CSC_idx      = [ 1:length( nNCsf ) ];
%     
%     %%
%     try
%         launch_par_preproc_EM(params);
%     catch ME
%         fid = fopen([params.savepath,'errorLog_',seshs(it).name,'.txt'],'w');
%         fprintf(fid,ME.message);
%         fclose(fid);
%     end;
%     
% end;
% 
% %% set params
% params                  = [];
% params.pID              = 'P08';
% params.p2Nlxdata = ['/media/rouxf/rds-share/Archive/MICRO/',params.pID,'/fVSpEM/'];%path 2 Nlx data
% params.savepath = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',params.pID,'/fVSpEM/'];%path for saving
% params.mode         = 'stimlocked';%
% params.pre          = 5;% baseline time range
% params.post         = 10;% params.post stim time range
% params.makeLFP      = 1;
% params.makeSpikes   = 1;
% params.makeMUA      = 0;
% params.rpath = '/home/rouxf/';
% params.Nlxdat = dir(params.p2Nlxdata);
% params.Nlxdat(1:2) = [];
% params.Nlxdat = {params.Nlxdat.name}';
% 
% %%
% seshs = dir([params.p2Nlxdata,'2017*']);
% {seshs(:).name}'
% 
% %% launch parallel preproc
% for it = 1:length(seshs)
%             
%     %%
%     params.sesh         = it;
%     nNCsf               = dir([params.p2Nlxdata,filesep,seshs(it).name,'/*.ncs']);
%     params.CSC_idx      = [ 1:length( nNCsf ) ];
%     
%     %%
%     try
%         launch_par_preproc_EM(params);
%     catch ME
%         fid = fopen([params.savepath,'errorLog_',seshs(it).name,'.txt'],'w');
%         fprintf(fid,ME.message);
%         fclose(fid);
%     end;
% end;

%% set params
params                  = [];
params.pID              = 'P09';
params.p2Nlxdata = ['/media/rouxf/rds-share/Archive/MICRO/',params.pID,'/fVSpEM/'];%path 2 Nlx data
params.savepath = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',params.pID,'/fVSpEM/'];%path for saving
params.mode         = 'stimlocked';%
params.pre          = 5;% baseline time range
params.post         = 10;% params.post stim time range
params.makeLFP      = 1;
params.makeSpikes   = 1;
params.makeMUA      = 0;
params.rpath = '/home/rouxf/';
params.Nlxdat = dir(params.p2Nlxdata);
params.Nlxdat(1:2) = [];
params.Nlxdat = {params.Nlxdat.name}';

%%
seshs = dir([params.p2Nlxdata,'2017*']);
{seshs(:).name}'

%% launch parallel preproc
for it = 3%1:length(seshs)
            
    %%
    params.sesh         = it;
    nNCsf               = dir([params.p2Nlxdata,filesep,seshs(it).name,'/*.ncs']);
    params.CSC_idx      = [ 1:length( nNCsf ) ];
    
    %%
    try
        launch_par_preproc_EM(params);
    catch ME
        fid = fopen([params.savepath,'errorLog_',seshs(it).name,'.txt'],'w');
        fprintf(fid,ME.message);
        fclose(fid);
    end;
end;

% %% set params
% params                  = [];
% params.pID              = 'P09';
% params.p2Nlxdata = ['/media/rouxf/rds-share/Archive/MICRO/',params.pID,'/cnEM/'];%path 2 Nlx data
% params.savepath = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',params.pID,'/cnEM/'];%path for saving
% params.mode         = 'stimlocked';%
% params.pre          = 5;% baseline time range
% params.post         = 10;% params.post stim time range
% params.makeLFP      = 1;
% params.makeSpikes   = 1;
% params.makeMUA      = 0;
% params.rpath = '/home/rouxf/';
% params.Nlxdat = dir(params.p2Nlxdata);
% params.Nlxdat(1:2) = [];
% params.Nlxdat = {params.Nlxdat.name}';
% 
% %%
% seshs = dir([params.p2Nlxdata,'2017*']);
% {seshs(:).name}'
% 
% %% launch parallel preproc
% for it = 1:length(seshs)
%             
%     %%
%     params.sesh         = it;
%     nNCsf               = dir([params.p2Nlxdata,filesep,seshs(it).name,'/*.ncs']);
%     params.CSC_idx      = [ 1:length( nNCsf ) ];
%     
%     %%
%     try
%         launch_par_preproc_EM(params);
%     catch ME
%         fid = fopen([params.savepath,'errorLog_',seshs(it).name,'.txt'],'w');
%         fprintf(fid,ME.message);
%         fclose(fid);
%     end;
% end;

%% close pool ad exit
delete(gcp);
exit;