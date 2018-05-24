%% set the path def
restoredefaultpath;
addpath(genpath('C:\toolbox\AnalysisFred\'));
addpath(genpath('C:\toolbox\Nlx\MatlabImportExport_v6.0.0\'));
addpath(genpath('C:\toolbox\osort-v3.0\osort-v3-rel\code\'));
addpath(genpath('C:\toolbox\chronux_2_12\'));
addpath('C:\toolbox\fieldtrip-20160309\');
ft_defaults;
%% open pool of workers
%parpool('SpmdEnabled',false);
clear;
clc;

%% set the data-paths directions
p2Nlxdata = 'C:\Data\Nlx\fVSp\';
p2logdat = 'C:\Experiments\EMpairs_v4_2016-1007\log\EM\';
p2logparams = 'C:\Experiments\EMpairs_v4_2016-1007\params\';
%% select session to analyze
ses = 1;
pre = 5;
post = 5;
chanSel = [2];
%% create the session labels of the Nlx data
Nlxdat = {};
Nlxdat{1} = 'P02_EM01_fVSp_2016-Jul-10_19-15-00\';
Nlxdat{2} = 'P02_EM02_fVSp_2016-Jul-11_13-43-29\';
Nlxdat{3} = 'P02_EM03_fVSp_2016-Jul-11_19-41-38\';
Nlxdat{4} = 'P02_EM05_fVSp_2016-Jul-13_14-15-12\';
Nlxdat{5} = 'P02_EM06_fVSp_2016-Jul-15_11-13-75\';

%% create the session labels of the log data
EMlogfile = {};
EMlogfile{1} = 'P02_fVSp_SE01_10-Jul-2016_18_45_5180100_LogFile_EMtask.txt';
EMlogfile{2} = 'P02_fVSp_SE02_11-Jul-2016_13_41_5562500_LogFile_EMtask.txt';
EMlogfile{3} = 'P02_fVSp_SE03_11-Jul-2016_19_40_3995000_LogFile_EMtask.txt';
EMlogfile{4} = 'P02_fVSp_SE05_13-Jul-2016_14_15_1228100_LogFile_EMtask.txt';
EMlogfile{5} = 'P02_fVSp_SE06_15-Jul-2016_11_13_754400_LogFile_EMtask.txt';

EMparams{1} = 'P02_fVSp_10-Jul-2016_19_8_5566000_params_completed.mat';
EMparams{2} = 'P02_fVSp_SE02_11-Jul-2016_14_3_4109400_params_completed.mat';
EMparams{3} = 'P02_fVSp_SE03_11-Jul-2016_20_3_4330400_params_completed.mat';
EMparams{4} = 'P02_fVSp_SE05_13-Jul-2016_14_36_3377300_params_completed.mat';
EMparams{5} = 'P02_fVSp_SE06_15-Jul-2016_11_37_1374700_params_completed.mat';

load([p2logparams,EMparams{ses}]);
%% get the Retrieval log data
params.p = p2logdat;
params.fn = EMlogfile{ses};
params.ntrl = 150;
params.ncol = 12;

[LogDat1] = getNewLogDataEM(params,'RET');

%% get indexes of conditions
k1 = 0;k2 = 0;k3 = 0;
pidx = [];
fidx = [];
fpidx = [];
for it = 1:size(LogDat1.log,1)
        
    cond ={LogDat1.log{it,3}(1) LogDat1.log{it,4}(1)};
    
    if strcmp([cond{:}],'pp')
        k1 = k1+1;
        pidx(k1) = it;
    elseif strcmp([cond{:}],'ff')
        k2 = k2+1;
        fidx(k2) = it;
    elseif strcmp([cond{:}],'fp')
        k3 = k3+1;
        fpidx(k3) = it;
    end;

end;

cond = sort([pidx fidx fpidx]);
%% get log data from Encoding
params.ncol = 9;
[LogDat2] = getNewLogDataEM(params,'ENC');

if size(LogDat1.log,1) ~= size(LogDat2.log,1)
    error('number of trials must match');
end;

enc_dur = str2double(LogDat2.log(:,end));

[n,x] = hist(enc_dur);
figure;bar(x,n);
%% get TTLs
filename = 'Events.nev';

[events,~] = getRawTTLs([p2Nlxdata,Nlxdat{ses},filename]);

[ttl_idx] = find(events(:,2)==7);

if ( length(LogDat1.log)+length(LogDat2.log) ~= length(ttl_idx) )
    error('number of ttls must match number of events');
end;

chck = [events(ttl_idx,2) events(ttl_idx,1) events(ttl_idx+1,2) events(ttl_idx+1,1)];
dt1 = diff(chck(:,[2 4]),[],2)./1e6;

dt2 = (events(ttl_idx+1,1)-events(ttl_idx,1))./1e6;

jitter1 = round((max(dt1)-min(dt1))*1000)/1000;
jitter2 = round((max(dt2)-min(dt2))*1000)/1000;

if (jitter1 > 0.005) || (jitter2 > 0.005) || (max(dt1) >2) || (max(dt2) >2) || ~isequal(jitter1,jitter2)
    error('stimulus on-off duration must remain below 2s and jitter below 5 ms');
end;
%% extract Fs and ADV
filename = 'CSC_LA1.ncs';

[timestamps,~,hdr] = getRawCSCData( [p2Nlxdata,Nlxdat{ses},filename], 1, 1, 1 );
FirstTimeStamp = timestamps(1);

[idx] = getNlxHeaderParam(hdr,'SamplingFrequency');
Fs = str2double(hdr{idx}(20:end));

[idx] = getNlxHeaderParam(hdr,'ADBitVolts');
ADV = str2double(hdr{idx}(12:end))*1e6;
%%
[hdr2] = ft_read_header([p2Nlxdata,Nlxdat{ses}]);
[event] = ft_read_event([p2Nlxdata,Nlxdat{ses}]);

[ttl_idx2] = find([event(:).value]==7);

if ( sum(diff([ttl_idx ttl_idx2'],[],2)) ~= 0 )
    error('wrong trigger aligment');
end;
%% convert ts to samples
timeStampsPerSample = 1/Fs*1e6;

samples = zeros(size(events,1),1);
samples2 = zeros(size(events,1),1);
parfor it = 1:size(events,1)
    timestamps;
    samples(it) = (events(it,1)-double(FirstTimeStamp))./timeStampsPerSample + 1;
    samples2(it) = (event(it).timestamp-double(hdr2.FirstTimeStamp))./hdr2.TimeStampPerSample + 1;
end;

sample3 = [event(ttl_idx2).sample]';

if unique(events(ttl_idx,2)) ~=7
    error('trigger values are out of range');
end;
%clear events;

%%
%samples = round(samples);
samples2 = samples2(ttl_idx);
eventts = [event(ttl_idx).timestamp]';
xt = str2double(LogDat2.log(:,5));

if length(samples2)/2 ~=length(LogDat1.log)
    error('number of events out of range');
end;
%% get encoding and retrieval trials
del_idx = [];
for it = 1:size(LogDat1.idx,1)
    LogDat1
    n = 0;
    if it > 1
        n = LogDat1.idx(it-1,2)*2;
    end;
    idx = n + 1: n+length(LogDat1.idx(it,1):LogDat1.idx(it,2));%ENC
    
    del_idx = [del_idx idx(end)+1:idx(end)+length(idx)];%RET
end;

%% sanity check
x = [];
for it = 1:size(LogDat1.idx,1)
    x = [x ones(1,length(LogDat1.idx(it,1):LogDat1.idx(it,2))) 2*ones(1,length(LogDat1.idx(it,1):LogDat1.idx(it,2)))];
end;
del_idx2 = find(x==2);

if ~isequal(del_idx,del_idx2)
    error('trial deletion indices out of range');
end;

%del_idx = find(x==1);

%% delete trials
samples2(del_idx) = [];
eventts(del_idx) = [];

if length(samples2) ~=length(LogDat1.log)
    error('number of events out of range');
end;
%% keep only succesfull trials
sel_idx = find(sum(str2double(LogDat1.log(:,5:6)),2)==2);
cond = cond(sel_idx);
samples2 = samples2(sel_idx);
eventts = eventts(sel_idx);
xt = xt(sel_idx);
params.iti1 = params.iti1(sel_idx);
enc_dur = enc_dur(sel_idx);

% if (sum(diff([diff(xt) diff(eventts)./1e6],[],2))) ~=0
%     error('event timing is out of sync');
% end;

%% segmentation indexes
trl = zeros(length(samples2),3);
for it =1:length(samples2)
    trl(it,:) = [samples2(it)-(pre*Fs) samples2(it)+(post*Fs) -(pre*Fs)];
end;

%% do spectral analysis for channel selection
S1 = cell(length(chanSel),1);
S2 = cell(length(chanSel),1);

files = dir([p2Nlxdata,Nlxdat{ses},'*.ncs']);

for ct = 1%1:1:length(chanSel)
        
    st = tic;
    fprintf([num2str(ct),'/',num2str(length(chanSel))]);
    fprintf('\n');
    
    filename = files(chanSel(ct)).name;
    
    [~,dataSamples,~] = getRawCSCData( [p2Nlxdata,Nlxdat{ses},filename], [], [], 1 );
    
    %[tsi] = timeStampinter(timestamps);
    
    dataSamples = dataSamples*ADV;

    %%
    datafilt  = dataSamples;
%     Wp = 500/Fs/2;
%     Ws = 1200/Fs/2;
%     [n,Wn] = buttord(Wp,Ws,3,60);
%     [b,a] = butter(n,Wn,'low');
%     
%     [datafilt] = filtfilt(b,a,dataSamples);    
%     clear a b;
%     clear dataSamples timestamps;
    
    %[datafilt] = locdetrend(datafilt,Fs,[.1 .05]);
    
    %[cleanSignal noise] = CleanLineNoise(datafilt','Fs',Fs,'noiseFreq',50,'windowSize',1);
    %clear datafilt;
    %clear noise;
    %%    
    x = zeros(size(trl,1),length(trl(1,1):trl(1,2)));
    parfor it = 1:size(trl,1)
        datafilt;
        trl;
        x(it,:) = datafilt(trl(it,1):trl(it,2));
    end;
    
    datafilt = x';%dims are time x trial
    clear x;
    
    %%
    p = max(abs(datafilt).^2,[],1);
    z = (p-mean(p))./std(p);
    datafilt(:,find(z>=2)) = [];
    params.iti1(find(z>=2)) = [];
    enc_dur(find(z>=2)) = [];
    %%
    [v,s_idx] = sort(params.iti1);
    [v2,s_idx2] = sort(enc_dur);
    yt = -v;
    yt2 =v2-7;
    
    figure;
    subplot(5,1,1:4);imagesc(-pre:1/32000:post,1:size(datafilt,2),datafilt(:,s_idx)');
    hold on;
    for jt = 1:length(yt);plot(yt(jt),jt,'k.');end;
    plot([0 0],[1 size(datafilt,1)],'Color',[.5 .5 .5]);
    subplot(515);plot(-pre:1/32000:post,mean(datafilt,2),'b-');
    
    figure;
    subplot(5,1,1:4);imagesc(-pre:1/32000:post,1:size(datafilt,2),datafilt(:,s_idx2)');
    hold on;
    for jt = 1:length(yt);plot(yt2(jt),jt,'k.');end;
    plot([0 0],[1 size(datafilt,1)],'Color',[.5 .5 .5]);
    subplot(515);plot(-pre:1/32000:post,mean(datafilt,2),'b-');
    %%
    ktapers = 8;
    NW = (ktapers+1)/2;
    
    params = [];
    params.tapers =  [NW ktapers];
    params.Fs = Fs;
    params.fpass = [0 Fs/2];
    params.pad = 4;
    
    datafilt2 = datafilt;    
    parfor it = 1:size(datafilt,2)        
        fprintf([num2str(it),'/',num2str(size(datafilt,2))]);
        [datafilt2(:,it)]=rmlinesc(datafilt2(:,it)',params,.05,'n');
        fprintf('\n');
    end;
    datafilt = datafilt2;
    clear datafilt2;  
    %%
    ktapers = 1;
    NW = (ktapers+1)/2;
    
    movingwin = [1 .25];
    
    params = [];
    params.Fs = Fs;
    params.pad = 8;
    params.tapers = [NW ktapers];
    params.fpass = [2 20];
    params.err = 0;
    params.trialave = 0;
    
    [S1{ct},t1,f1]=mtspecgramc(datafilt,movingwin,params);
    
    %%
    ktapers = 8;
    NW = (ktapers+1)/2;
    
    movingwin = [.25 .067];
    
    params = [];
    params.Fs = Fs;
    params.pad = 4;
    params.tapers = [NW ktapers];
    params.fpass = [20 100];
    params.err = 0;
    params.trialave = 0;
    
    [S2{ct},t2,f2]=mtspecgramc(datafilt,movingwin,params);
    
    clear datafilt;    

    toc(st);
end;

%%    
clear samples;
clear ttl_idx;
clear trl;

%%
figure;
for it = 1:length(chanSel)
    subplot(round(length(chanSel)/8)+1,8,it);
    imagesc(t1-pre,f1,squeeze(mean(20*log10(S1{it}),3))');axis xy;
    title(files(chanSel(it)).name);
end;

%%
figure;
for it = 1:length(chanSel)
    subplot(round(length(chanSel)/8)+1,8,it);
    imagesc(t2-post,f2,squeeze(mean(20*log10(S2{it}),3))');axis xy;
    title(files(chanSel(it)).name);
end;