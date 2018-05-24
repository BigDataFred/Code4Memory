function preprocEM( pID , sesh )

%% set the session parameters
if nargin ==0
    pID = 'P04';
    sesh = 1:3;
end;

%% get 32 CPUs to do parallel processing of channels
if isempty(gcp('nocreate'))
    parpool(32,'SpmdEnabled',false);
end;

%% set path env
restoredefaultpath;
addpath(genpath('~rouxf/releaseDec2015/'));%needed to read Nlx data
addpath(genpath('~rouxf/osort-v3-rel/'));%needed for spike detection & sorting
addpath(genpath('~rouxf/AnalysisFred/'));
addpath(genpath('~rouxf/wave_clus-testing/'));
addpath('~rouxf/fieldtrip-20161009/');
ft_defaults;

%%
for pt = 1:numel(sesh)
    
    %% get the session data info
    [sdat] = create_session_dat(pID,sesh(pt));
    
    %% set baseline and poststim duration (always add 1xtime the tf window length)
    pre = 3;%baseline    
    post = 3;%post stim

    %% get the logfile data
    params1.p = sdat.p2lf;%path
    params1.fn = sdat.lfn;%filename
    params1.ntrl = 49;%number of trials
    params1.ncol = 9;% number of columns in logfile
    
    %same as above
    params2.p = sdat.p2lf;
    params2.fn = sdat.lfn;
    params2.ntrl = 49;
    params2.ncol = 12;
    
    %read logfile for Encoding
    [LogDat1,trl_idx1] = getNewLogDataEM(params1, 'ENC');
    %read logfile for Retrieval
    [LogDat2,~] = getNewLogDataEM(params2, 'RET');
    
    [RTs] = str2double(LogDat1.log(:,end));%get RTs
    post = post+RTs;
    
    %find hits, misses and stim categories
    cat = LogDat1.log(:,3:4);
    c = {};
    parfor jt = 1:size(cat,1)
        
        c{jt} = [cat{jt,1}(1) cat{jt,2}(1)];% reads out info from logdata
        
    end;
    
    ix    = {};
    ix{1} = find(strcmp(c,'fp'));% indices corresponding to face-place tr
    ix{2} = find(strcmp(c,'pp'));% indices corresponding to place-place tr
    ix{3} = find(strcmp(c,'ff'));% indices corresponding to face-face tr
    
    ix{4} = find(sum([str2double(LogDat2.log(:,5:6))],2)==2);% both images were correctly remebered
    ix{5} = find(sum([str2double(LogDat2.log(:,5:6))],2)==1);% only 1/2 images were correctly remebered
    ix{6} = find(sum([str2double(LogDat2.log(:,5:6))],2)==0);%no images were correctly remembered
    
    ix_readme = {'face-place' 'place-place' 'face-face' 'correct-both' 'miss-one' 'miss-both'};% save info for later
    LogDat_readme = {'LogDat1:Enc','LogDat2:Ret'};
    
    %%    save the logfile data
    if isempty(dir([sdat.p2Nlxf,filesep,'log_dat',filesep]))
        mkdir([sdat.p2Nlxf,filesep,'log_dat',filesep]);
    end;
    [~,fn,~] = fileparts(sdat.lfn);
    
    save([sdat.p2Nlxf,filesep,'log_dat',filesep,fn,'_LogDat.mat'],'LogDat1','LogDat2','ix','trl_idx1','ix_readme','LogDat_readme','RTs');
        
    %% check for existing dir
    if isempty(dir([sdat.p2Nlxf,filesep,'spike_dat',filesep]))
        mkdir([sdat.p2Nlxf,filesep,'spike_dat',filesep]);
    end;
    
    %% check for existing dir
    if isempty(dir([sdat.p2Nlxf,filesep,'mua_dat',filesep]))
        mkdir([sdat.p2Nlxf,filesep,'mua_dat',filesep]);
    end;
    
    %% check for existing dir
    if isempty(dir([ sdat.p2Nlxf,filesep,'lfp_dat',filesep ]))
        mkdir([ sdat.p2Nlxf,filesep,'lfp_dat',filesep ]);
    end;
    
    %% read in the CSC header
    files = dir([sdat.p2Nlxf,'*CSC_*.ncs']);% get the filenames
    
    [timestamps,~,nrSamples,sampleFreq,~,headerInfo] = getRawCSCTimestamps( [sdat.p2Nlxf,files(1).name] );% read in header-info
    
    hdr = [];
    hdr.Fs = sampleFreq;
    
    [idx] = getNlxHeaderParam(headerInfo,'ADBitVolts');% conversion factor
    hdr.scalef = str2double(headerInfo{idx}(12:end));
    
    [idx] = getNlxHeaderParam(headerInfo,'AcqEntName');% CSC name
    hdr.label = headerInfo{idx}(13:end);
    hdr.nChans = length(hdr.label);
    
    hdr.FirstTimeStamp = int64(timestamps(1));% first ts
    hdr.nSamples = nrSamples;
    hdr.TimeStampPerSample = (timestamps(2)-timestamps(1))/512;
    
    %% read the TTLs
    [event] = ft_read_event([sdat.p2Nlxf,'Events.nev']);
    
    for it = 1:length([event(:).value])
        event(it).samples = double((event(it).timestamp - uint64(hdr.FirstTimeStamp))/uint64(hdr.TimeStampPerSample+ 1));
    end;
    
    ttl_idx = find([event(:).value] == 7);
    ttl_idx = ttl_idx([trl_idx1{:}]);
    
    if length(ttl_idx) ~= length(LogDat1.log)
        error('trial and trigger mismatch');
    end;
    
    parfor kt = 1:length(event)
        event(kt).samples = double((event(kt).timestamp - uint64(hdr.FirstTimeStamp))/uint64(hdr.TimeStampPerSample+ 1) );
    end;
    
    s = [event(ttl_idx).samples]';
    
    %% make the trl matrix
    %[TRL] = int64([s-pre*hdr.Fs s+post*hdr.Fs -pre*hdr.Fs*ones(length(s),1)]);
    [TRL] = int64([s-pre*hdr.Fs s+post*hdr.Fs -pre*hdr.Fs*ones(length(s),1)]);
    %% pre-allocate
    lfp_data    = cell(1,length(files));
    spike_data1  = cell(1,length(files));
    spike_data2  = cell(1,length(files));
    spike_data3  = cell(1,length(files));
    spike_data4  = cell(1,length(files));
    mua_data    = cell(1,length(files));
    
    %% spike-band filter settings
    Wp = [300 3000]./(hdr.Fs/2);
    [b,a] = butter(4,Wp);
    
    Hd{1} = b;
    Hd{2} = a;

    %% segment the data
    parfor it = 1:length(files)% parloop over channels
        
        %% tmp variables
        hdr2 = hdr;
        sdat2 = sdat;
        timestamps2 = timestamps;
        files2 = files;
        
        [~,chan,~] = fileparts( files2(it).name );

        %% read the CSC data
        % indexes for reading in the .ncs data
        fromInd = timestamps2(1);
        toInd = timestamps2(end);                
        
        [~,dataSamples] = getRawCSCData( [sdat2.p2Nlxf,files2(it).name], fromInd, toInd , 4);
        dataSamples = double(dataSamples(:))';
        dataSamples = dataSamples.*hdr2.scalef.*1e6;
        
        %% estimate threshold parameters for clustering
        [par]                               = set_parameters_Bham(hdr2.Fs);
        par.fnamespc                        = [ par.fnamespc,chan ] ;
        par.fname_in                        = [ par.fname_in,chan ];
        
        [fs,~,~,~] = amp_detect(dataSamples,par);
        
        waveclus                            = [];
        waveclus.thr                        = par.sort_fmin * median(abs(fs))/0.6745;
        waveclus.thrmax                     = par.stdmax * median(abs(fs))/0.6745;
        fs = [];
        
        %% make the spike times
        nsegs = 6;
        
        lts = length(timestamps2);
        segl = floor(lts/nsegs);
        
        tsmin = 1: segl : lts;
        tsmin = tsmin(1: nsegs );
        tsmax = tsmin - 1;
        tsmax = tsmax(2:end);
        tsmax = [tsmax, lts];
        
        tsmin = timestamps2(int64(tsmin));
        tsmax = timestamps2(int64(tsmax));
        
        %% pre-allocate
        ts1 = [];        thr1 = [];        wavf1 = [];  
        ts2 = [];        thr2 = [];        wavf2 = []; noiseTraces = [];    wavfo = [];     blockOffsets = [];  stdEstimates =[];
        for zt = 1:length(tsmin)%% spike detection
            
            fromInd = tsmin(zt);
            toInd = tsmax(zt);
            [ts_raw,dataSamples2] = getRawCSCData( [sdat.p2Nlxf,files2(it).name], fromInd, toInd , 4);
            
            x = double(dataSamples2(:))';
            x = x*hdr.scalef*1e6;
            
            %% spike detection using wave_clus                                    
            [~,spikeWaveforms1,thr1,spikeTimestamps1] = amp_detect(x,par);
            spikeTimestamps1 = spikeTimestamps1*1e6/sampleFreq+tsmin(zt);%convert samples to timestamps and add start time of system clock
            
            close(gcf);
            
            ts1 = [ts1 spikeTimestamps1];
            wavf1 = [wavf1;spikeWaveforms1];
            thr1 = [thr1 thr1];
            
            %% spike detection using osort
            params_ts = [];
            params_ts.nrNoiseTraces = 0;
            params_ts.detectionMethod= 1;% 1 -> from power signal, 2 threshold positive, 3 threshold negative, 4 threshold abs, 5 wavelet
            %params.detectionParams: depends on detectionMethod.
            if params_ts.detectionMethod==1; params_ts.kernelSize = 18;end;
            if params_ts.detectionMethod==4; params_ts.scaleRanges = [];params_ts.waveletName = '';end;
            params_ts.alignMethod = 1;% 1-> find peak, 2->none, 3->peak of power signal, 4->peak of MTEO signal.
            params_ts.extractionThreshold = 5;
            params_ts.prewhiten = 0;
            params_ts.samplingFreq = sampleFreq;
            params_ts.limit = 1000;
            
            [rawMean, filteredSignal, rawTraceSpikes,spikeWaveforms2, spikeTimestamps2, runStd2, upperlim, noiseTracesTmp] = extractSpikes(x', Hd, params_ts );
            
            wavf2o = spikeWaveforms2;%original waveform
            
            stdEstimates(zt) = std(filteredSignal);
            
            if size(noiseTracesTmp,1)>1
                noiseTracesTmp2=[];
                noiseTracesTmp2(1:size(noiseTracesTmp,1) ,1) = ones( size(noiseTracesTmp,1), 1 )*i;
                noiseTracesTmp2(1:size(noiseTracesTmp,1),2:size(noiseTracesTmp,2)+1)=noiseTracesTmp;
                noiseTraces = [noiseTraces; noiseTracesTmp2];
            end
            
            %upsample and re-align
            spikeWaveforms2=upsampleSpikes(spikeWaveforms2);
            spikeWaveforms2 = realigneSpikes(spikeWaveforms2, spikeTimestamps2, params_ts.alignMethod, stdEstimates(zt));  %3==type is negative, do not remove any if too much shift
            
            %convert timestamps
            if ~isempty(spikeTimestamps2)
                spikeTimestamps2 = convertTimestamps( ts_raw, spikeTimestamps2, params_ts.samplingFreq, 2 );
            end;
            
            blockOffsets = [blockOffsets ts_raw(1)];
                        
            wavf2 = [wavf2;spikeWaveforms2];
            wavfo = [wavfo;wavf2o];
            ts2 = [ts2 spikeTimestamps2];
            thr2 = [ params_ts.extractionThreshold ];
            
            
        end;
        
        %whiten
        if ( size(noiseTraces,1)>0 ) && ( size(wavfo,1)>0 )
            [~, transUp, corr, ~] = posthocWhiten(noiseTraces, wavfo, params_ts.alignMethod);
            wavfo = transUp; transUp = [];
            
            %only store the autocorrelation,not all noise traces
            noiseTraces=corr; corr = [];
        end;
        
        %% 
        
        waveclus.spikes                     = wavf1;  
        waveclus.index                      = ts1; 
        
        [sortedSpikes1] = doSpikeSorting_waveclus( waveclus , par );
        waveclus = [];
        
        %%
        osort                               = [];
        osort.ADbitVolts                    = hdr2.scalef;
        osort.allSpikes                     = wavf2;    wavf2 = [];
        osort.allSpikesCorrFree             = wavfo;    
        osort.allSpikesTimestamps           = ts2;     
        osort.blocksProcessed               = zt;
        osort.noiseTraces                   = noiseTraces;
        osort.blockOffsets                  = blockOffsets;
        osort.dataSamplesRawUncorrected     = [];
        osort.dataSamplesRaw                = dataSamples2;
        osort.stdEstimates                  = stdEstimates;
        osort.upperlim                      = upperlim;
        osort.runStd2                       = runStd2;
        osort.rawTraceSpikes                = rawTraceSpikes;
        osort.filteredSignal                = filteredSignal;
        osort.rawMean                       = rawMean;
        
        [sortedSpikes2] = doSpikeSorting_osort( osort );
        osort = [];
        rawMean = []; filteredSignal= []; rawTraceSpikes= [];spikeWaveforms2= []; spikeTimestamps2= []; runStd2= []; upperlim= []; noiseTracesTmp= [];
        %% make dummy ft structure for sorted spike data
        wft = linspace(1/sampleFreq,64*1/sampleFreq,64);
                
        %timestamps defined with waveclus                   
        [spike1]                     =  struct;
        spike1.label                 = {chan};
        spike1.timestamp{1}          = int64(sortedSpikes1.newSpikeTimes);
        spike1.waveform{1}(1,:,:)    = sortedSpikes1.wavf(:,1:2:end)';
        spike1.unit{1}               = sortedSpikes1.assignedCluster;
        spike1.hdr                   = hdr;
        spike1.dimord                = '{chan}_lead_time_spike';
        spike1.cfg.thr               = thr1;
        spike1.waveformtime          = wft(1:2:end);
        
        sortedSpikes1 = [];
        
        %timestamps defined with O-sort
        [spike2]                     =  struct;
        spike2.label                = {chan};        
        spike2.timestamp{1}          = int64(sortedSpikes2.newSpikesTimestampsNegative);
        spike2.waveform{1}(1,:,:)    = sortedSpikes2.newSpikesNegative(:,1:8:end)';
        spike2.unit{1}               = sortedSpikes2.assignedClusterNegative;        
        spike2.hdr                   = hdr;
        spike2.dimord                = '{chan}_lead_time_spike';
        spike2.cfg.thr               = thr2;
        spike2.waveformtime          = wft(1:2:end);
        
        sortedSpikes2 = [];
        
        %% make dummy ft structure for unsorted spike data
        
        %timestamps defined with waveclus
        [spike3]                     =  struct;
        spike3.label                 = {chan};        
        spike3.timestamp{1}          = int64(ts1);
        spike3.waveform{1}(1,:,:)    = wavf1(:,1:2:end)';
        spike3.unit{1}               = zeros(1,length(ts1));
        spike3.hdr                   = hdr;
        spike3.dimord                = '{chan}_lead_time_spike';
        spike3.cfg.thr               = thr1;
        spike3.waveformtime          = wft(1:2:end);
        
        thr1 = [];
        ts1 = [];
        wavf1 = [];
        
        %timestamps defined with O-sort
        [spike4]                     =  struct;
        spike4.label                 = {chan}; 
        spike4.timestamp{1}          = int64(ts2);
        spike4.waveform{1}(1,:,:)    = wavfo(:,1:2:end)';
        spike4.unit{1}               = zeros(1,length(ts2));
        spike4.hdr                   = hdr;
        spike4.dimord                = '{chan}_lead_time_spike';
        spike4.cfg.thr               = thr2;
        spike4.waveformtime          = wft(1:2:end);
        
        thr2 = [];
        ts2 = [];
        wavfo = [];
        %% do the spike time segmentation
        
        cfg                         = [];
        cfg.trl                     = [ TRL ];
        cfg.trlunit                 = 'samples';
        cfg.hdr                     = hdr;
        
        [spike_data1{it}]            = ft_spike_maketrials(cfg , spike1);
        [spike_data2{it}]            = ft_spike_maketrials(cfg , spike2);
        spike1 = [];
        spike2 = [];
        
        cfg                         = [];
        cfg.trl                     = [ TRL ];
        cfg.trlunit                 = 'samples';
        cfg.hdr                     = hdr;
        
        [spike_data3{it}]            = ft_spike_maketrials(cfg , spike3);
        [spike_data4{it}]            = ft_spike_maketrials(cfg , spike4);
        spike3 = [];
        spike4 = [];
        %% create dummy structure with CSC data        
        
        % interpolate the timestamps to make time axis continuous
        [tsi] = timeStampinter(timestamps2);%interpolate timestamps
        tsi = tsi-min(tsi);
        tsi = tsi./1e6;
        
        %dummy ft structure
        CSC = [];
        CSC.fsample = hdr2.Fs;
        CSC.label = {chan};
        CSC.trial{1} = dataSamples;
        CSC.time{1} = tsi;
        
        tsi = [];
        dataSamples = [];
        
%         %% make the MUA                
%         
%         %bandpass spike-band
%         cfg                     = []; 
%         cfg.bpfilter            = 'yes';
%         cfg.bpfreq              = [500 5000];
%         cfg.bpfiltord           = 2;
%         cfg.padding             = 20;
%         cfg.padtype             = 'mirror';
%                 
%         [dum] = ft_preprocessing(cfg , CSC );        
%         
%         %take the envelope
%         cfg                     = [];
%         cfg.hilbert             = 'abs';
%         
%         [dum] = ft_preprocessing(cfg, dum );
%         
%         % lowpass at 200 Hz
%         cfg                     = [];
%         cfg.lpfilter            = 'yes';
%         cfg.lpfreq              = 200;
%         cfg.lpfilttype          = 'but';
%         cfg.lpfiltord           = 2;
%         
%         [dum] = ft_preprocessing(cfg , dum );
%         
%         % segment trials
%         cfg                     = [];
%         cfg.trl                 = [ TRL ];
%         
%         [dum] = ft_redefinetrial(cfg, dum );
%         
%         %downsample
%         cfg                     = [];
%         cfg.resamplefs          = 1000;
%         cfg.detrend             = 'no';
%         
%         [dum] = ft_resampledata(cfg , dum ); 
%         
%         % notch filter line noise
%         bsf = [ 50:50:200 ];
%         for kt = 1:length(bsf)
%             cfg                     = [];
%             cfg.bsfilter            = 'yes';
%             cfg.bsfreq              = [bsf(kt)-2 bsf(kt)+2];
%             
%             [dum] = ft_preprocessing(cfg , dum );            
%         end;
%                 
%         [mua_data{it}] = dum;
%         dum = [];
%         
%         %% preprocessing of LFP data
%         cfg                     = [];
%         cfg.padding             = 20;
%         cfg.padtype             = 'mirror';
%         cfg.lpfilter            = ['yes'];
%         cfg.lpfreq              =  300;
%         
%         [CSC] = ft_preprocessing(cfg, CSC );
%         
%         %% do the LFP segmentation
%         
%         cfg                 = [];
%         cfg.trl             = [ TRL ];
%         
%         [CSC] = ft_redefinetrial(cfg, CSC );
%         
%         %%
%         cfg                         = [];
%         cfg.detrend                 = 'no';
%         cfg.resamplefs              = 1000;
%         
%         [CSC]                  = ft_resampledata(cfg , CSC );        
%         
%         %%
%         signal = CSC.trial;
%         CSC.trial{1} = {};
%         for zt = 1:length(CSC.label)
%             for yt = 1:length(signal)
%                 [cleanSignal,~] = CleanLineNoise( signal{yt}(zt,:) ,'Fs', CSC.fsample , 'noiseFreq', 50,'windowSize',1);
%                 signal{yt}(zt,:) = cleanSignal;
%             end;
%             cleanSignal = [];
%         end;
%         CSC.trial = signal;
%         signal = [];
%         
%         %%
%         lfp_data{it} = CSC;
        CSC = [];

    end;
    
    %% concatenate spike data across channels
    
    dum = {};
    [dum{1}] = ft_appendspike([],spike_data1{:});
    [dum{2}] = ft_appendspike([],spike_data2{:});
    
    [spike_data] = dum;
    clear spike_data1 spike_data2 dum;
    
    %% save spike data to disc
    spk_read_me = {'waveclus','osort'};
    
    save([sdat.p2Nlxf,filesep,'spike_dat',filesep,fn,'_sortedSPIKEdata.mat'],'spike_data','spk_read_me');
    clear spike_data;
    
    %% concatenate spike data across channels
    
    dum = {};
    [dum{1}] = ft_appendspike([],spike_data3{:});
    [dum{2}] = ft_appendspike([],spike_data4{:});
    
    [spike_data] = dum;
    clear spike_data3 spike_data4 dum;
    
    %% save spike data to disc
    spk_read_me = {'waveclus','osort'};
    
    save([sdat.p2Nlxf,filesep,'spike_dat',filesep,fn,'_unsortedSPIKEdata.mat'],'spike_data','spk_read_me');
    clear spike_data;
    
    %% concatenate spike data across channels    
    [mua_data] = ft_appenddata([],mua_data{:});
    
    %% save spike data to disc
    mua_readme = {'mua:Self'};
    
    save([sdat.p2Nlxf,filesep,'mua_dat',filesep,fn,'_preprocMUAdata.mat'],'mua_data','mua_readme');
    clear mua_data*;
    
    %% concatenate LFP data across all channels    
    [lfp_data] = ft_appenddata([],lfp_data{:});% concatenate BF data across al channels
     
    n = [];
    for jt = 1:length(lfp_data.time);
        n(jt) = max(lfp_data.time{jt});
    end;
    %% save LFP-data to disc
    lfp_readme = {'Fs:1kHz','50Hz-noise free'};
    save([sdat.p2Nlxf,filesep,'lfp_dat',filesep,fn,'_preprocLFPdata_ds.mat'],'lfp_data','lfp_readme');
    
    clear lfp_data;
    
end;
delete(gcp);%shutdown parpool
exit;%exit MATLAB
