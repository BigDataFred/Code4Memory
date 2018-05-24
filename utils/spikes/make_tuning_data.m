function make_tuning_data(make_figure,chan_number,sf,logf)
if nargin ==0
    sf.p2d = '/media/rouxf/My Passport/';
    sf.p2sd = [sf.p2d,'A/sort/5/'];
    make_figure = 'yes';
    chan_number = [ 46 ];
    %logf = [];
    logf = '/home/rouxf/Data/EM/Log/P01/LogMat/P01_25-May-2016_10_13_3612100_log_tune.mat';
end;
%% set path env
restoredefaultpath;
addpath('~rouxf/MATLAB/toolboxes/fieldtrip-20160309/');
ft_defaults;
addpath(genpath('~rouxf/MATLAB/utils/'));
addpath(genpath('~rouxf/MATLAB/toolboxes/releaseDec2015/'));
addpath(genpath('~rouxf/MATLAB/toolboxes/osort-v3-rel/'));
%% load the sorted data
spdat = load([sf.p2sd,'A',num2str(chan_number),'_sorted_new.mat']);
%% read Nlx header and get TTLs
try
    [hdr]   = ft_read_header(sf.p2d);
    [event] = ft_read_event(sf.p2d);
catch
    filename = [sf.p2d,'A/A',num2str(chan_number),'.Ncs'];
    [hdr] = extract_header_info( filename );
end;

chck1 = double(((double(hdr.LastTimeStamp-hdr.FirstTimeStamp))/1e6)/60);
chck2 = double(hdr.nSamples/hdr.Fs/60);

if ( round(chck2)-round(chck1) ~= 0 )
    error('timepoints and sample number must match');
end;
%% adjust time stamps from absolute system time to recording time
st = double(spdat.newTimestampsNegative-double(hdr.FirstTimeStamp))/1e6;
%% spike train autocorrelation
AC = cell(1,length(spdat.useNegative));
for it = 1:length(spdat.useNegative)
    sel_idx = find(spdat.assignedNegative == spdat.useNegative(it));
    [AC{it}]= compute_spiketrain_autocorr(st(sel_idx));
end;
%% use for sanity check
%st2 = [et et+0.01 et+0.1 et+0.2 et+0.3];
%st = sort(st2);
%% return trigger indices
cID = {0};
if ~isempty(logf)
    load(logf);
    ntrl = length(LogDat.dat(2:end,1));
    [ttl_idx] = check_ttl_vs_log(LogDat.params,[event(:).value],ntrl);
    et = double([event(:).timestamp]-double(hdr.FirstTimeStamp))/1e6;% adjust time stamps from absolute system time to recording time
    cID = unique(LogDat.cond);
end;
%%
dataset = [sf.p2d,'A/A',num2str(chan_number),'.Ncs'];
fn = dir(dataset);
if isempty(fn)
    fprintf('Warning:file not detected\n');
    return;
end;
%% segment the data into trigger-aligned trials
if cID{1} ~=0
    spike_dat = cell(length(spdat.useNegative),length(cID));
    
    params.Fs = hdr.Fs;
    params.bw = 0.001;
    params.pre = -1;
    params.post = 1;
    params.window = 'hanning';
    
    for it = 1:length(spdat.useNegative)
        
        sel_idx = find(spdat.assignedNegative == spdat.useNegative(it));
        
        for jt = 1:length(cID)
            
            sel_idx2 = find(ismember(ttl_idx,find(strcmp(LogDat.cond,cID(jt)))));
            
            [spike_dat{it,jt}] = segment_spike_data(et(sel_idx2),st(sel_idx),params);
            
        end;
    end;
    
    
%     %% trigger based LFP segmentation
%     
%     x = double([event(ttl_idx).sample]);
%     x = x';
%     
%     [trl]             = [x-1*hdr.Fs x+1*hdr.Fs -1*hdr.Fs*ones(size(x,1),1)];
%     del_idx = [find(trl(:,2) > hdr.nSamples) find(trl(:,1) < 1)];
%     trl(del_idx,:) = [];
%     
%     [CSC_preproc] = launch_LFP_preproc(dataset,trl);
%     
%     %compute spectrum of trigger aligned LFP
%     cfg = [];
%     cfg.method = 'mtmfft';
%     cfg.pad = 'maxperlen';
%     cfg.foi = 1:300;
%     cfg.taper = 'dpss';
%     cfg.tapsmofrq = 1;
%     
%     [fft] = ft_freqanalysis(cfg,CSC_preproc);
%     
%     %compute event-related LFP
%     
%     cfg = [];
%     cfg.preproc.demean = 'yes';
%     cfg.preproc.lpfilter = 'yes';
%     cfg.preproc.lpfreq = 30;
%     cfg.preproc.lpfilttype = 'but';
%     cfg.preproc.padtype = 'mirror';
%     cfg.preproc.padding = CSC_preproc.time{1}(end)*5;
%     
%     [tlck] = ft_timelockanalysis(cfg,CSC_preproc);
%     clear CSC_preproc;
end;
% %% spike based LFP segmentation
% fft2  = cell(1,length(spdat.useNegative));
% tlck2 = cell(1,length(spdat.useNegative));
% for it = 1:length(spdat.useNegative)
%     
%     [sel_idx] = find(spdat.assignedNegative == spdat.useNegative(it));
%     
%     x = zeros(length(sel_idx),1);
%     for jt = 1:length(sel_idx)
%         x(jt) = ( spdat.newTimestampsNegative(sel_idx(jt)) - double(hdr.FirstTimeStamp) )./hdr.TimeStampPerSample + 1;
%     end;
%     
%     [trl]             = [x-1*hdr.Fs x+1*hdr.Fs -1*hdr.Fs*ones(size(x,1),1)];
%     del_idx = [find(trl(:,2) > hdr.nSamples) find(trl(:,1) < 1)];
%     trl(del_idx,:) = [];
%     
%     [CSC_preproc2] = launch_LFP_preproc(dataset,trl);
%     
%     %compute spectrum of trigger aligned LFP
%     
%     cfg = [];
%     cfg.method = 'mtmfft';
%     cfg.pad = 'maxperlen';
%     cfg.foi = 1:6000;
%     cfg.taper = 'dpss';
%     cfg.tapsmofrq = 1;
%     
%     %[fft2{it}] = ft_freqanalysis(cfg,CSC_preproc2);
%     
%     %compute event-related LFP
%     
%     cfg = [];
%     cfg.preproc.demean = 'yes';
%     %cfg.preproc.lpfilter = 'yes';
%     %cfg.preproc.lpfreq = 30;
%     %cfg.preproc.lpfilttype = 'but';
%     %cfg.preproc.padtype = 'mirror';
%     %cfg.preproc.padding = CSC_preproc2.time{1}(end)*5;
%     
%     [tlck2{it}] = ft_timelockanalysis(cfg,CSC_preproc2);
%     clear CSC_preproc2;
%     
% end;
%% visualize results
if strcmp(make_figure,'yes')
    %%
    visualize_sorting(spdat);
    %%
    C = {};
    C{1} = [.9 0 0];
    C{2} = [0 0 .9];
    C{3} = [0 0 0];
    
    figure;
    a1 = zeros(size(spike_dat,1),1);
    a2 = zeros(size(spike_dat,1),1);
    for it = 1:size(spike_dat,1)
        for jt = 1:size(spike_dat,2)
            
            subplot(2,3,it);
            a1(it) = gca;
            hold on;
            rasterPlot(spike_dat{it,jt}.trlst,C{jt});
            
            subplot(2,3,it+3);
            a2(it) = gca;
            hold on;
            visualize_pSTH(spike_dat{it,jt}.dt,spike_dat{it,jt}.frate,C{jt});
            
        end;
    end;
    set(a1,'Xlim',([min(spike_dat{1,1}.dt) max(spike_dat{1,1}.dt)]));
    
    for it =1:length(a1)
        title(a1(it),num2str(spdat.useNegative(it)));
    end;
    
    yl = zeros(length(a1),1);
    for it =1:length(a2)
        yl(it) = max(get(a2(it),'YTick'));
    end;
    set(a2,'YLim',[0 max(yl)]);
    %%
    figure;
    k = 0;
    for it = 1:length(AC)
        k= k+1;
        subplot(2,ceil(length(AC)/2),k);
        visualize_spikeAC(AC{it});
        ylim([0 100]);
    end;
    %%
    figure;
    subplot(221);
    hold on;
    plot([min(tlck.time) max(tlck.time)],[0 0],'r--');
    plot(tlck.time,tlck.avg);
    subplot(222);
    hold on;
    plot([min(fft.freq) max(fft.freq)],[0 0],'r--');
    plot(fft.freq,20*log10(fft.powspctrm));
    subplot(223);
    hold on;
    plot([min(tlck2{1}.time) max(tlck2{1}.time)],[0 0],'r--');
    plot(tlck2{1}.time,tlck2{1}.avg);
    subplot(224);
    hold on;
    plot([min(fft2{1}.freq) max(fft2{1}.freq)],[0 0],'r--');
    plot(fft2{1}.freq,20*log10(fft2{1}.powspctrm));
    %%
    for it = 1:get(gcf,'Number')
        set(it,'Color','w');
    end;
end;

%%
return















