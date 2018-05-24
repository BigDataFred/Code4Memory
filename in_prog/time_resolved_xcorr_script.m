%% set path defs
restoredefaultpath;
addpath('/media/rouxf/rds-share/Fred/code/mcode/custom/project_EM/');
addpath('/media/rouxf/rds-share/Common/fieldtrip-20170115/');
ft_defaults;
addpath(genpath('/media/rouxf/rds-share/Common/chronux_2_11/'));

%% data parameters 
pID = 'P07';% subject IDcd 
Exp = 'fVSpEM';% exp mode (can be fVSpEM or cnEM)
seshSel =1;% number of the session to run analysis on

%% get the labels of each session for selected data-set
rp = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',pID,'/',Exp,'/'];% root-path

sesh = dir(rp);% search for directories
sesh(1:2) = [];% eliminate empty info

% select only those directories that contain preprocessed NLX data
x = {};
for it = 1:length( sesh )
    x(it) = {sesh(it).name};
end;
sesh = x;
chck = regexp(sesh,'\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}');%searches for pattern 
chck = [chck{:}];% binarize dir name
sesh(chck==0)=[];% retain only those dirs that match pattern
sesh{seshSel}%output to command line

%% get files with spk-data & lfp-data
p2spkd          = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',pID,'/',Exp,'/',sesh{seshSel},'/spike_dat/'];
[spkFiles]      = dir([p2spkd,'*.mat'])%spk

p2lfp           = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',pID,'/',Exp,'/',sesh{seshSel},'/lfp_dat/'];
[lfpFiles]      = dir([p2lfp,'*.mat'])%lfp

%% load the log-data
p2log           = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',pID,'/',Exp,'/',sesh{seshSel},'/log_dat/'];
[logf]          = dir([p2log,'*_LogFile_EMtask_LogDat.mat'])%log
load([p2log,logf.name])

%% extract trials corresponding to encoding phase
[trl_all]       = extract_Enc_trl(LogDat1)%output to command line

%% init some vars for later usage
trlID           = trl_all; % trial labels for spike times
nCSC            = length( spkFiles);% n of channels
nTrl            = length(trlID);% n of trials
lag             = 250;% lag up to which xcorr will be computed

[XC]            = zeros(nCSC,nCSC,lag*2+1);% matrix with xcorr for channel pairs
[spkDat]        = cell(1,nCSC);% cell array with spike data, ft-format
[lfpDat]        = cell(1,nCSC);% cell array with lfp data, ft-format
[ts]            = cell(1,nCSC);% cell array with spike times
[spk]           = cell(1,nCSC);% cell array with spike times

%% init cell arrays to store lfp and spike times for each channel
[data1]         = cell(1,nCSC); % samples x trials matrix for lfp data
[data2]         = cell(1,nCSC); % samples x trials matrix for spike times

%% time bins for the spike times
dta             = -2e3:1:5e3;% whole trial epoch
dtb             =  2e3:1:5e3;% encoding time window

%% load the spk-data
for it = 1:nCSC% 1st loop over CSC channels
    
    fprintf([num2str(it),'/',num2str(nCSC)]);%monitor progress        
    
    %% load the spike data
    dat             = load([p2spkd,spkFiles(it).name]);
    [spkDat{it}]    = dat.save_data{1}{1}{3};% extract spike data in ft-format
       
    %% extract and convert the spike times
    ts{it} = spkDat{it}.time{1}.*1e3;% convert spike times to ms
    
%     %% load the lfp data
%     dat = load([p2lfp,lfpFiles(it).name]);
%     [lfpDat{it}] = dat.save_data{1}{1}{1};
%     
%     %%
%     cfg                     = [];
%     cfg.trials              = trl_all;
%     cfg.latency             = [dta(1) dta(end)];
%     
%     [lfpDat{it}] = ft_selectdata( cfg , lfpDat{it} );
%     
%     %%
%     cfg                     = [];
%     cfg.keeptrials          = 'yes';
%     
%     [tlck] = ft_timelockanalysis( cfg , lfpDat{it} );
%     
%     %%
%     data1{it} = squeeze(tlck.trial);
%     data1{it} = data1{it}';    
%     clear tlck;
    
    %%
    fprintf('\n');
    
end;

%% organize spike times into binary samples x trials matrix
for it = 1:nCSC% 1st loop over CSC channels
    
    fprintf([num2str(it),'/',num2str(nCSC)]);%monitor progress
    
    %% loop over trials and bin spike times in ms time bins
    dum1 = zeros( length(trl_all) , length(dta) );
    dum2 = zeros( length(trl_all) , length(dtb) );
    
    x1 = spkDat{it}.trial{1};
    x2 = ts{it};
    for kt = 1:nTrl
                
        dum3 = x2(x1 == trlID(kt));% readout of spike times for trial of interest
        dum3(dum3<dta(1) & dum3>dta(end)) = [];
        
        [dum1(kt,:)] = histc(dum3,dta);
        
        dum3 = x2(x1 == trlID(kt));% readout of spike times for trial of interest
        dum3(dum3<dtb(1) & dum3>dtb(end)) = [];
        
        [dum2(kt,:)] = histc(dum3,dtb);
        
    end;
    dum1(:,end) = [];%get rid of last sample
    data2{it} = dum1';% store spike times for later usage
    
    dum2(:,end) = [];%get rid of last sample
    spk{it} = dum2';% store spike times for later usage
    
    clear x1 x2 dum1 dum2 dum3;
    %%
    fprintf('\n');
    
end;

%% compute pairwise x-correlation
nc = zeros(nCSC,nCSC);
for it = 1:nCSC% 1st loop over CSC channels
    
    fprintf([num2str(it),'/',num2str(nCSC)]);%monitor progress                 
    
    x1 = spk{it};
    
    %% loop over channels
    xc = zeros(nCSC,lag*2+1);
    for lt = 1:nCSC% 2nd loop over channels
        
        x2 = spk{lt};    cnt = 0;
        
        for kt = 1:nTrl
            
            [spk1] = x1(:,kt);% extract the binned spike times
            [spk2] = x2(:,kt);% extract the binned spike times
            
            fr1 = sum(spk1)/3;    fr2 = sum(spk2)/3;% compute firing rates
            
            if (fr1 >5) %&& (fr2 >5)
                cnt =  cnt +1;    
                xc(lt,:) = xc(lt,:) + xcorr(spk1',spk2',lag);
            end;
            
        end;
        xc(lt,:) = xc(lt,:)./cnt;    
        nc(it,lt) = cnt;
    end;
    XC(it,:,:) = xc;
    
    %%
    fprintf('\n');
    
end;

%%
for it = 1:size(XC,1);
    c = 0;
    figure;
    for jt= 1:size(XC,1);
        c=c+1;
        subplot(nCSC/8,8,c);
        hold on;
        h =bar(-lag:lag,squeeze(XC(it,jt,:)));
        set(h,'FaceColor','k');
        plot([-lag lag],[2 2],'r--');
        axis tight;
        ylim([0 15]);
        title([nc(it,jt)]);
    end,
end;

%%
movingwin               = [1 0.025];

T = movingwin(1);
W = 5;
TW = T*W;
k = 2*TW-1;

params                  = [];
params.Fs               = 1e3;%lfpDat.fsample;
params.fpass            = [30 150];
params.pad              = 2;
params.err              = [1 0.05];
params.trialave         = 1;
params.tapers           = [TW k];

[Pnn,t,f,Serr]=mtspecgramc(data1{2},movingwin,params);

ix = find(t-2 < 0);

[Pxx,t,f,R,Serr]=mtspecgrampb(data2{5},movingwin,params,0);
R = ones(length(f),1)*R(:,1)';


[C,phi,S12,S1,S2,t,f]=cohgrampb(data2{2},data2{5},movingwin,params,0);

figure;
subplot(211);
hold on;
imagesc(t-2,f,abs(S1)');axis xy;
plot([0 0],[min(f) max(f)],'w');
plot([2 2],[min(f) max(f)],'w');
axis tight;

B  = mean(Pnn(ix,:));
SD = std(Pnn(ix,:),0,1);
B  = ones(length(t),1)*B;
SD = ones(length(t),1)*SD;
 
%Z = (Pnn-B)./SD;

figure;
subplot(211);
hold on;
imagesc(t-2,f,(Pxx'./R));axis xy;
plot([0 0],[min(f) max(f)],'w');
plot([2 2],[min(f) max(f)],'w');
axis tight;
subplot(212);
hold on;
imagesc(t-2,f,20*log10(Pnn)');axis xy;
plot([0 0],[min(f) max(f)],'w');
plot([2 2],[min(f) max(f)],'w');
axis tight;

%%
XC2 = {};    raster = {};    cnt = 0;

for it = [13]%1:length( spkFiles)% loop over CSC channels       
    
    %% loop over the clusters
    for zt = 1:length( spkDat{it}.unit)
         
        cnt =cnt+1;                 
        
        %% loop over trials
        xc = zeros(lag*2+1,length(dta));
        cnt2 = 0;
        cnt3 = 0;
        for jt = 1:length(trlID)
            
            fprintf([num2str(jt),'/',num2str(length(trlID))]);
            
            cnt2 = cnt2+1;            
            ix = find(spkDat{it}.trial{zt} == trlID(jt));% find spike times corresponding to trial
            
            if ( ~isempty(ix) ) && ( length(ix) > 1 )
                
                ts2 = ts{it}(ix);% readout of spike times for trial of interest
                 
                [spk3] = histc(ts2,dta);% do the binning of the spike times
                
                selIx = find(dta>=2e3 & dta <=5e3);
                fr = sum(spk3(selIx))/(length(selIx)/1e3);
                
                raster{cnt}{cnt2} = ts2;
                
                if fr >5     
                    cnt3 = cnt3+1;
                    for kt = 1:length( dta ) % loop over bins                        
                        if (kt>lag) && (kt < length( spk3 )-lag)
                            ix2 = kt-lag/2:kt+lag/2;% shift window for autocorrelation in 1 ms steps
                            dum = xcorr(spk3(ix2));% compute autocorrelationx
                            xc(:,kt) = xc(:,kt) +dum';% sum autocorr across trials for each bin
                        end;                        
                    end;                    
                end;
                
            end;
            fprintf('\n');
            
        end;
        xc =  xc./cnt3;% normalize by the number of trials
                
        %if max(max(xc)) >=.1
        XC2{cnt} = xc;% if autocorrelation exceeds threshold save matrix
        %end;
    end;
    
end;

%% visualize the time-resolved autocorrelation

dtc = -2000:100:5000;
for it = 1:length( XC2 )
    
    figure;
    subplot(5,2,[1 2 3 4]);
    hold on;
    n = [];
    for jt = 1:length(raster{it})
        x = raster{it}{jt};
        if ~isempty(x)
            n(jt,:) = histc(x,dtc);
            y = jt*ones(1,length(x));
            x = [x;x];
            y =[y+.25;y-.25] ;
            h = line(x,y);
            set(h,'Color','k');
        end;
    end;
    plot([0 0],[0 jt+1],'r');
    plot([2e3 2e3],[0 jt+1],'r');
    ylim([0 jt+1]);
    
    subplot(5,2,[5 6]);
    bar(dtc,sum(n,1)/size(n,1)/0.1,'k','LineWidth',3);
    axis tight;
    
    subplot(5,2,[7 8 9 10]);
    hold on;
    imagesc(dta,-lag:lag,XC2{it});
    plot([0 0],[-lag lag],'w');
    plot([2e3 2e3],[-lag lag],'w');
    %caxis([0 0.1]);
    axis tight;
    
end;
