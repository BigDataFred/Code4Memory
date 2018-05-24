function [p2SpkDat,spkFiles,p2micDat,micFiles,p2macDat,macFile,macDat,ixIEDdat,logDat] = setupMicANDMacANDLogPathANDFileNames4EMtaskAMS(pID,expMode,expSel,sesh,seshSel,BFlabel,BFlabel2,bfSel,bfSel2,MWlabel,mwSel,micDatMode)

%% load the Logfile
p2logDat = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',pID,'/',expMode{expSel},'/',sesh{seshSel},'/log_dat/'];
x = dir(p2logDat);
chck = [];
for it = 1:length(x)
    chck(it) = ~isempty(regexpi(x(it).name,'LogFile_EMtask_LogDat.mat'));
end;

logFile = dir([p2logDat,x(chck==1).name]);

[logDat] = load([p2logDat,logFile.name]);


%% search for the spike data
p2SpkDat = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',pID,'/',expMode{expSel},'/',sesh{seshSel},'/spike_dat/'];
spkDat = [pID,'_TS',num2str(seshSel),'_spike_data_*',BFlabel2{bfSel2},MWlabel{mwSel},'_',sesh{seshSel},'_stimlocked.mat'];
[spkFiles] = dir([p2SpkDat,spkDat]);

%% searching the MW-LFP-data
p2micDat = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',pID,'/',expMode{expSel},'/',sesh{seshSel},'/lfp_dat/'];
if micDatMode == 2
    micDat = [pID,'_',expMode{expSel},'_*',BFlabel{bfSel},'_',sesh{seshSel},'_preprocANDrejectedIED_stimlocked.mat'];
else
    micDat = [pID,'_TS',num2str(seshSel),'_lfp_data_*',BFlabel{bfSel},'*_',sesh{seshSel},'_stimlocked.mat'];
end;
[micFiles] = dir([p2micDat,micDat]);

% %% search for Macro-LFP-data
p2macDat = [];
macFile = [];
macDat = [];
% p2macDat = ['/media/rouxf/rds-share/iEEG_DATA/MACRO/',pID,'/',expMode{expSel},'/',sesh{seshSel},'/'];
% macFile = dir([p2macDat,pID,'_PreprocMacroLFP_fVSpEM_',sesh{seshSel},'.mat']);
% if ~isempty( macFile )
%     [macDat] = load([p2macDat,macFile.name]);
% else
%     [macDat] = [];
% end;

%% load IED-based trial rejection indexes
ixIEDdat = [];
ixIEDdat.LFPdatMW = {};
ixIEDdat.delIx = [];
ixIEDdat.samp = [];
ixIEDdat.selIx = [];
IEDfile = dir([p2micDat,pID,'_',expMode{expSel},'_*',BFlabel{bfSel},'_',sesh{seshSel},'_preprocANDrejectedIED_stimlocked.mat']);
if ~isempty(IEDfile)
    IEDfile = dir([p2micDat,pID,'_',expMode{expSel},'_*',BFlabel{bfSel},'_',sesh{seshSel},'_preprocANDrejectedIED_stimlocked.mat']);
elseif isempty(bfSel)
    for it  = 1:length(BFlabel)
        IEDfile = dir([p2micDat,pID,'_',expMode{expSel},'_*',BFlabel{it},'_',sesh{seshSel},'_preprocANDrejectedIED_stimlocked.mat']);
        if ~isempty(IEDfile)
            [dum] = load([p2micDat,IEDfile.name]);
            [ixIEDdat.LFPdatMW{it}] =[];
            [ixIEDdat.delIx] = [ixIEDdat.delIx dum.delIx];
            [ixIEDdat.samp] = [];
            [ixIEDdat.selIx] = [];
        end;
    end;
    ixIEDdat.delIx = unique(ixIEDdat.delIx);  
else
    ixIEDdat = [];
end;

%%
%%
if isempty(ixIEDdat.delIx)

    IEDfile = dir([p2micDat,pID,'_',expMode{expSel},'_',sesh{seshSel},'_IEDtrlIdx4Rejection_stimlocked.mat']);
    [dum] = load([p2micDat,IEDfile.name]);
    [ixIEDdat.LFPdatMW] =[];
    [ixIEDdat.samp] = [];
    [ixIEDdat.selIx] = [];
    
    %dum.sel(2) = {1:50};% for testing
    %dum.sel(4) = {1:50};% for testing
    
    delIx = dum.sel;
    for it = 1:length(delIx)
        if size(delIx{it},1) > size(delIx{it},2);delIx{it} = delIx{it}';end;
    end;
    
    n = [];for it = 1:length(delIx);n(it) = length(delIx{it});end;
    [~,sIdx] = sort(n);
    
%     x = [];for it = 1:length(delIx); x = [x delIx{it}];end;
%     x(x==0) = [];
    
    pct1 = length([delIx{:}])/length(logDat.LogDat1.log);
    if pct1 < 0.5
        [ixIEDdat.delIx] = unique([delIx{:}]);
        [ixIEDdat.chanDel] = [];
        [ixIEDdat.chanName] = {};
        ixIEDdat.delIx(ixIEDdat.delIx==0)=[];
    else
        pass = [];
        nopass = [];
        x = [];cnt1 = 0; cnt2 = 0;
        for it = 1:length(sIdx)
            pct2 = length(unique([x delIx{sIdx(it)}]))/length(logDat.LogDat1.log);
            if pct2 < 0.5
                x = [x delIx{sIdx(it)}];
                cnt1 = cnt1+1;
                pass(cnt1) = sIdx(it);
            else
                cnt2 = cnt2+1;
                nopass(cnt2) = sIdx(it);
            end;
        end;
        if ~isempty(pass)
            [ixIEDdat.delIx] = unique([delIx{pass}]);
            [ixIEDdat.chanDel] = unique([delIx{nopass}]);
            ixIEDdat.delIx(ixIEDdat.delIx==0)=[];
            ixIEDdat.chanDel(ixIEDdat.chanDel==0)=[];
            ixIEDdat.chanName = dum.chanName(nopass);
        else
            [ixIEDdat.delIx] = unique([delIx{:}]);
            [ixIEDdat.chanDel] = [];
            ixIEDdat.delIx(ixIEDdat.delIx==0)=[];
            ixIEDdat.chanName = [];
        end;
    end;     
    
    if ~isempty(ixIEDdat.chanName)
        sel = [];cnt = 0;
        for jt = 1:length(ixIEDdat.chanName)
            for it = 1:length(spkFiles)
                if ~isempty(cell2mat(regexpi(spkFiles(it).name,ixIEDdat.chanName(jt))));
                    cnt = cnt+1;
                    sel(cnt) = it;
                end;
            end;
        end;
        spkFiles(sel) = [];
        
        sel = [];cnt = 0;
        for jt = 1:length(ixIEDdat.chanName)
            for it = 1:length(micFiles)
                if ~isempty(cell2mat(regexpi(micFiles(it).name,ixIEDdat.chanName(jt))));
                    cnt = cnt+1;
                    sel(cnt) = it;
                end;
            end;
        end;
        micFiles(sel) = [];
        
    end;
    
end;

