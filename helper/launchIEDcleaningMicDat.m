%% set the path environment
addpath('/home/rouxf/tbx/fieldtrip-20170618/');
ft_defaults;

addpath('/home/rouxf/prj/Bham/code/mcode/project_EM/');
addpath('/home/rouxf/prj/Bham/code/mcode/utils/LFP/');
addpath('/home/rouxf/prj/Bham/code/mcode/visualize/continuous/');

%% open parallel pool
if isempty(gcp('nocreate'))
    parpool(12,'SpmdEnabled',false);
end;

%% get the patient IDs
rpath = '/media/rouxf/rds-share/iEEG_DATA/MICRO/';
x = dir(rpath);

pID = cell(1,length(x));
for it = 1:length( x )
    pID(it) = {x(it).name};
end;

pID(1:2) = []

%% loop over the patients
for pt = 1:length( pID )
    
    %% get the exp-modes
    ppath = [rpath,pID{pt},filesep];
    x = dir(ppath);
    x(1:2) = [];
    expMode = cell( 1 , length(x) );
    for it = 1:length(x)
        expMode(it) = {x(it).name};
    end;
    
    chck1 = regexpi(expMode,'[fVSpEM]');
    chck2 = regexpi(expMode,'[cnEM]');
    n1 = zeros(length(chck1),1);
    n2 = zeros(length(chck1),1);
    for it = 1:length(chck1)
        n1(it) = length( chck1{it} );
        n2(it) = length( chck2{it} );
    end;

    selIx = unique([find(n1==6);find(n2==4)]);
    expMode = expMode(selIx)
    
    %% loop over the experiment modes
    for et = 1:length(expMode)
        
        %% get the session labels
        spath = [ppath,expMode{et},filesep];
        x = dir(spath);
        x(1:2) = [];
        sesh = cell(1,length(x));
        for it = 1:length(x)
            sesh(it) = {x(it).name};
        end;
        sesh
        
        %% loop over the different sessions
        for st = 1:length( sesh )
            
            %%
            [p2MIC] = [spath,sesh{st},filesep,'lfp_dat',filesep];
            
            [micFiles] = dir([p2MIC,'*lfp_data*stimlocked.mat']);
            
            if ~isempty( micFiles )
                
                [segMode] = micFiles(1).name(max(regexp(micFiles(1).name,'_'))+1:end-4);
                
                %% load the MW data
                [micDat] = cell( 1 , length(micFiles) );
                parfor it = 1:length(micFiles)
                    
                    fprintf([num2str(it),'/',num2str( length(micFiles) )]);
                    [dat] = load([p2MIC,micFiles(it).name]);
                    [micDat{it}] = dat.save_data{1}{1}{1};
                    fprintf('\n');
                    
                end;
                
                [micDat] = ft_appenddata( [] , micDat{:} );
                
%                 %% preproc Micro-data
                cfg                     = [];
                cfg.detrend             = 'yes';
                cfg.demean              = 'yes';
                cfg.bpfilter            = 'yes';
                cfg.bpfilttype          = 'but';
                cfg.bpfiltord           = 2; 
                cfg.bpfreq              = [20 50];
%                 cfg.padtype             = 'mirror';
%                 cfg.padding             = 5;
                %cfg.hilbert             = 'abs';
                
                [dum] = ft_preprocessing( cfg , micDat );
                
                %% average the signals of the MWs of each BF-shank
                [dum] = makeAVG4BF(dum);
                
                selIx = [];
                for it = 1:length(dum.label)
                    selIx(it) = ~isempty(dum.label{it});
                end;
                
                cfg                     =[];
                cfg.channel             = dum.label(selIx==1);
                
                [dum] = ft_selectdata( cfg ,dum );
                
                cfg                     = [];
                cfg.latency             = [-1 5];
                
                [dum] = ft_selectdata(cfg, dum);
                
                %%
                nsamp = length(dum.trial{1});
                trl = [];
                for it = 1:length(dum.trial)
                    trl(it,:) = [(it-1)*nsamp+1 (it-1)*nsamp+nsamp 0];
                end;
                
                sel = {};artifact = {};
                for ct = 1:length(dum.label)
                    cfg                     = [];
                    cfg.trl                 =  trl;
                    cfg.continuous          = 'no';
                    cfg.artfctdef.zvalue.channel = dum.label(ct);
                    cfg.artfctdef.zvalue.cutoff = 4;
                    cfg.artfctdef.zvalue.interactive = 'no';
                    [cfg, artifact{ct}] =  ft_artifact_zvalue(cfg, dum);
                                                            
                    sel{ct}(1) = 0;
                    cnt = 0;
                    for ft = 1:size(trl,1)
                        for ht = 1:size(artifact{ct},1)
                            d = diff([artifact{ct}(ht,1) artifact{ct}(ht,2)])/1e3;
                            if d >= 0.015 && d <= 0.075
                                if any(ismember([artifact{ct}(ht,1) artifact{ct}(ht,2)],trl(ft,1):trl(ft,2)))
                                    cnt = cnt+1;
                                    sel{ct}(cnt) =ft;
                                end;
                            end;
                        end;
                    end;
                    if ~isempty(sel{ct})
                        sel{ct} = unique(sel{ct});
                    end;
                end;
                
                %                 %% preproc Micro-data
                cfg                     = [];
                cfg.detrend             = 'yes';
                cfg.demean              = 'yes';
                cfg.bpfilter            = 'yes';
                cfg.bpfilttype          = 'but';
                cfg.bpfiltord           = 2; 
                cfg.bpfreq              = [1 35];
%                 cfg.padtype             = 'mirror';
%                 cfg.padding             = 5;
                %cfg.hilbert             = 'abs';
                
                [dum2] = ft_preprocessing( cfg , micDat );
%                 
                                %% average the signals of the MWs of each BF-shank
                [dum2] = makeAVG4BF(dum2);
                
                selIx = [];
                for it = 1:length(dum2.label)
                    selIx(it) = ~isempty(dum2.label{it});
                end;
                
                cfg                     =[];
                cfg.channel             = dum2.label(selIx==1);
                
                [dum2] = ft_selectdata( cfg ,dum2 );
                
                cfg                     = [];
                cfg.latency             = [-1 5];
                
                [dum2] = ft_selectdata(cfg, dum2);
                
                %% detect the IED events
                [chanIx,trlIx,ixON,ixOFF] = artifactRejectIED(dum,dum2);
                
                %%
                [trl,pct] = IEDdetector(dum,dum.label,[-1 5]);
                for it = 1:length(trl);if size(trl{it},2)>size(trl{it},1);trl{it} = trl{it}';end;end;
                for it = 1:length(sel);if size(sel{it},2)>size(sel{it},1);sel{it} = sel{it}';end;end;    
                
                %%
                for it = 1:length(sel);sel{it} = sel{it}(ismember(sel{it},trlIx{it}));end;
                for it = 1:length(sel);if isempty(sel{it});sel{it}=0;end;end;
                for it = 1:length(sel);sel{it} = unique([sel{it};trl{it}]);end;
                                
                savename = [pID{pt},'_',expMode{et},'_',sesh{st},'_IEDtrlIdx4Rejection_',segMode,'.mat'];
                chanName = dum.label;
                save([p2MIC,savename],'sel','chanName');
                
%                 %% smaple indexes across trials
%                 nsamp = length(micDat.trial{1});
%                 ntrl = length( micDat.trial );
%                 
%                 samp = zeros(ntrl,2);
%                 for it = 1:ntrl
%                     samp(it,1) = (it-1)*nsamp+1;
%                     samp(it,2) = (it-1)*nsamp+nsamp;
%                 end;
%                 
%                 %% detect the IED events
%                 [chanIx,trlIx,ixON,ixOFF] = artifactRejectIED(dum,dum2);
%                 
%                 %% reject trials with IEDs for each channel individually
%                 for it = 1:length( micDat.label )
%                     
%                     delIx = trlIx{it};
%                     selIx = setdiff(1:length(micDat.trial),delIx);
%                     
%                     cfg                     = [];
%                     cfg.trials              = selIx;
%                     cfg.channel             = micDat.label(it);
%                     
%                     [LFPdatMW] = ft_selectdata( cfg , micDat );
%                     
%                     savename = [pID{pt},'_',expMode{et},'_',micDat.label{it},'_',sesh{st},'_preprocANDrejectedIED_',segMode,'.mat'];
%                     save([p2MIC,savename],'LFPdatMW','selIx','delIx','samp');
%                     
%                 end;
                
            end;

            
        end;
    end;
end;
delete(gcp);
exit;