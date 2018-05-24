function preproc_macro_data_EM(pID,expMode)

%% default
if nargin ==0
    pID = 'P04';
    expMode = 'fVSpEM';
end;

%% set the path
restoredefaultpath;
ft = dir(['/home/rouxf/tbx/fieldtrip-********']);
addpath(['/home/rouxf/tbx/',ft.name]);
ft_defaults;
addpath('/home/rouxf/prj/Bham/code/bash/');
addpath('/home/rouxf/prj/Bham/code/mcode/project_EM/');
addpath('/home/rouxf/prj/Bham/code/mcode/utils/LFP/');

%% find different sessions
r_path = ['/media/rouxf/rds-share/Archive/MICRO/',pID,'/',expMode,'/'];

[sesh] = extract_sesh_label(r_path);

%% do the signal processing
path2MAC     = ['/media/rouxf/rds-share/Archive/MACRO/',pID,'/EDF/'];
savepath4MAC = ['/media/rouxf/rds-share/iEEG_DATA/MACRO/',pID,'/',expMode,'/'];
if ~isdir(savepath4MAC);
    [s] = mkdir(savepath4MAC);
end;

%%
cnt = 0;% just a counter
for it = length(sesh)%1:length(sesh)%loop over the different sessions        
    
    [micTS] = sesh(it).name;% name of microwire data
    [macFN] = find_macro_filename_EM(micTS,pID);% search lookup table for corresponding macro data
    
    if ~isempty(macFN)
        try
            
            cnt = cnt +1;% increment each time a macro-dataset is found
            
            %% check if time range of macro and micro data matches
            [hdr] = ft_read_header([path2MAC,macFN]);
            
            str = cell(1,length( hdr.orig.T0 ));
            for gt = 1:length(hdr.orig.T0)
                str(gt) = {sprintf('%02d',hdr.orig.T0(gt))};
            end;
            
            [macTS] = [str{1},'-',str{2},'-',str{3},'_',str{4},'-',str{5},'-',str{6}];
            
            if ~strcmp(micTS(1:10),macTS(1:10))
                error('date of micro and macro timestamps is out of range');
            end;            
            if ~strcmp(micTS(1:10),macTS(1:10))
                error('date of micro and macro timestamps is out of range');
            end;
            
            at = micTS(12:end);at(regexp(at,'-')) = ':';
            bt = macTS(12:end);bt(regexp(bt,'-')) = ':';                        
            
            et = etime(datevec(at),datevec(bt));
            if ( sign(et) == -1 )
                error('temporal offset between micro and macro data must be  positive');
            end;
            
             %% load the logfile information
            path2LOG = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',pID,'/',expMode,'/',sesh(it).name,'/log_dat/'];
            logFN = dir([path2LOG,'*LogDat*.mat']);
            load([path2LOG,logFN.name]);
            nTRL = size(LogDat1.log,1) +  size(LogDat2.log,1);% total number of trials in that session
            
            [trl_Enc] = extract_Enc_trl(LogDat1);
            [trl_Ret] = setdiff(1:nTRL,trl_Enc);
            
            rmbrd = ix{4}; % save the trials that were later rmbrd during retrieval
            frgtn = [ix{5};ix{6}];% save the trials that were later frgtn during retrieval
                                    
            % get the reaction times
            RTsRet = str2double(LogDat2.log(:,7));
            RTsEnc = str2double(LogDat1.log(:,end));
            
            [~,sIdx1] = sort(RTsEnc);% sort trials according to RT speeds
            [~,sIdx2] = sort(RTsRet);
            
            [dum1] = str2double(LogDat1.log(:,6))-str2double(LogDat1.log(:,5));
            [dum2] = str2double(LogDat2.log(:,11))-str2double(LogDat2.log(:,10));
            
            stimDT2 = [];
            stimDT2(trl_Enc) = dum1;
            stimDT2(trl_Ret) = dum2;
            
            %% read the trigger information
            cfg                     = [];
            cfg.dataset             = [path2MAC,macFN];
            cfg.channel             = {'TRIG'};
            
            [trigDat] = ft_preprocessing(cfg);% read the trigger channel from file            
            
            %% detect the trigger events
            dt = [0 diff(trigDat.trial{1})];% difference between adjacent samples
            [ixON] = find(sign(dt)==-1);% onsets of trigger
            [ixOFF] = find(sign(dt)==1);% offsets of trigger           
            
            [trigSampE,trigSampR,recT] = reconstructTTLtimingFromLog(ixON,ixOFF,trigDat,LogDat1,LogDat2);
            trigSamp = [trigSampE trigSampR]';
            delIx = find(sign(fliplr(recT'))==-1);
            
            %[stimDt] = diff([ixON' ixOFF'],[],2)/trigDat.fsample;
            %trlIdx = find( stimDt >1.95);% find trigger events that correspond to stimuli
            %delIdx = find(diff(trlIdx)~=1)+1;% these triggers are bogus
            %trlIdx(delIdx) = [];% delete bogus triggers            
            
            
%             trlVal = trigDat.trial{1}([trigSamp]);
%             trlVal = (trlVal-mean(trlVal))./std(trlVal);
%             trsh = quantile(abs(trlVal),0.95);
%             delIdx = find(abs(trlVal) >= trsh);
%             
%             if ~isempty(delIdx)
%                 trigSamp(delIdx)  = [];
%             end;
%             stimDt = stimDt(trlIdx);
            
           
            %% quick sanity check, always good to check your sanity ... it's
            % free and doesn't cost much ;)
            if ( length(trigSamp) ~= nTRL )
                error('WARNING: number of TTL and LOG events are not equal\n');
            end;
            
            % measure some delays
            presT = str2double(LogDat1.log(:,5:7));
            dt1 = presT(:,2)-presT(:,1); % delay between cue and pair
            dt2 = presT(:,3)-presT(:,2); % delay between pair and response
            
            % find samples for trigger events
            onSamp =trigSamp;            
            preSamp = 2*trigDat.fsample;% samples to substract to get the baseline
            postSamp = 5*trigDat.fsample;% samples to add to get the end of the poststim window
            
            % segment the macro data into epochs
            cfg                     =[];
            cfg.continuous          = 'no';
            cfg.dataset             = [path2MAC,macFN];
            cfg.channel             = {'*H*' , '-D*', '-Event', '-TRIG', '-OSAT', '-PR', '-Pleth' };
            cfg.bpfilter            = 'yes';
            cfg.bpfilttype          = 'fir';
            cfg.bpfreq              = [0.5 200];
            cfg.padtype             = 'data';
            cfg.padding             = 10;
            cfg.trl                 = round([onSamp-preSamp onSamp+postSamp -preSamp*ones(length(onSamp),1)]);
            cfg.demean              = 'yes';
            cfg.detrend             = 'yes';
            
            delIx = unique([find( cfg.trl(:,1) <=0) find(cfg.trl(:,2) > length(trigDat.trial{1}))]);
            cfg.trl(delIx,:) = [];
            
            [macroLFPdat] = ft_preprocessing( cfg );                        
            
            x = [];for zt = 1:length(macroLFPdat.time);x(zt,:) = [min(macroLFPdat.time{zt}) max(macroLFPdat.time{zt})];end;
            
            if min(min(round(x))) ~= -2
                error('trial time is out of range');
            end;
            
            if max(max(round(x))) ~= 5
                error('trial time is out of range');
            end;
            
            trlIdx_Enc = trl_Enc;
            trlIdx_Ret = trl_Ret;
            
            readme = {'macroLFPdat: basic preproc, no ref',...
                'trlIdx_Enc: indexes for Encoding',...
                'trlIdx_Ret: indexes for Retrieval',...
                'delIx: indexes of deleted trials',...
                'rmbrd: indexes for later remembered trials',...
                'frgtn: indexes for later forgotten trials'};
            
            %% save data to disk
            if ~isdir([savepath4MAC,sesh(it).name,filesep]);
                [s] = mkdir([savepath4MAC,sesh(it).name,filesep]);
            end;
            
            [savename] = [pID,'_PreprocMacroLFP_',expMode,'_',sesh(it).name,'.mat'];
            
            save([savepath4MAC,sesh(it).name,filesep,savename],'macroLFPdat','trlIdx_Enc','trlIdx_Ret','delIx','rmbrd','frgtn','readme');
            
        catch ME            
            fid = fopen([savepath4MAC,pID,'_PreprocMacroLFP_',expMode,'_',sesh(it).name,'_errorLog_b.txt'],'w+');
            fprintf(fid,ME.message);
            it
            sesh(it).name
            fprintf(ME.message);
        end;
    else
        fid = fopen([savepath4MAC,pID,'_PreprocMacroLFP_',expMode,'_',sesh(it).name,'_errorLog_b.txt'],'w+');
        fprintf(fid,['macro data for timestamp ',sesh(it).name,' is missing.']);
    end;
    
end;

%%
exit;