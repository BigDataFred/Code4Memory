function [data,readme] = analyze_EM_data(params)

data = {}; readme = {};

%%
try
        
    %% read the CSC data
    FieldSelection(1) = 1;%timestamps
    FieldSelection(2) = 0;
    FieldSelection(3) = 0;%sample freq
    FieldSelection(4) = 0;
    FieldSelection(5) = 1;%samples
    ExtractHeader = 0;
    
    ExtractMode = 1; % 2 = extract record index range; 4 = extract timestamps range.
    
    [params.timestamps, dataSamples] = Nlx2MatCSC_v3([params.p2Nlxdata,params.CSClabels{params.chanSel}], FieldSelection, ExtractHeader, ExtractMode, []);
    %flatten
    dataSamples=double(dataSamples(:))';
    dataSamples = dataSamples.*params.hdr.scalef.*1e6;
       
    %%
    if params.interrupted
        dx = diff( params.timestamps );
        ix = find(dx == max(dx))+1;
        [params.timestamps, dataSamples] = Nlx2MatCSC_v3([params.p2Nlxdata,params.CSClabels{params.chanSel}], FieldSelection, ExtractHeader, ExtractMode, []);
        params.timestamps = params.timestamps(ix+1:end);
        dataSamples = dataSamples(:,ix+1:end);
        if length(dataSamples) ~= length(params.timestamps)
            error('number of timestamps and samples must match');
        end;
        dataSamples=double(dataSamples(:))';
        dataSamples = dataSamples.*params.hdr.scalef.*1e6;
    end;
    
    %% tmp variables
    [~,chan,~] = fileparts( params.CSClabels{params.chanSel} );
    params.hdr2 = params.hdr;
    timestamps2 = params.timestamps;
    
    %% estimate threshold parameters for clustering
    [par]                               = set_parameters_Bham(params.hdr2.Fs);
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
    
    %% make the spike data
    if params.makeSpikes ==1
        
        %% pre-allocate
        ts1 = [];        thr1 = [];        wavf1 = []; noise_std =[]; spkIdx = {}; nsamp ={};
        ts2 = [];        thr2 = [];        wavf2 = []; noiseTraces = [];    wavfo = [];     blockOffsets = [];  stdEstimates =[];
        for zt = 1:length(tsmin)%% spike detection
            
            fromInd = tsmin(zt);
            toInd   = tsmax(zt);
            
            FieldSelection(1) = 1;%timestamps
            FieldSelection(2) = 0;
            FieldSelection(3) = 0;%sample freq
            FieldSelection(4) = 0;
            FieldSelection(5) = 1;%samples
            ExtractHeader = 0;
            
            ExtractMode = 4; % 2 = extract record index range; 4 = extract timestamps range.
            
            [ ts_raw, dataSamples2 ] = Nlx2MatCSC_v3([params.p2Nlxdata,params.CSClabels{params.chanSel}], FieldSelection, ExtractHeader, ExtractMode, [fromInd toInd]);
         
            x = double(dataSamples2(:))';
            x = x*params.hdr.scalef*1e6;
            
            %% spike detection using wave_clus
            [~,spikeWaveforms1,thr1,spikeTimestamps1,noise_std_detect,noise_std_sorted] = amp_detect(x,par);
            %spikeTimestamps1 = spikeTimestamps1.*1e6./params.hdr.Fs+tsmin(zt);%convert samples to timestamps and add start time of system clock                       
            
            spkIdx{zt} = spikeTimestamps1;% spkIdx+(nsamp)];% this will server to retrieve spikes times in continuous signal for LPF-spike interpolation
            nsamp{zt} = length(x);
            
             if ~isempty(spikeTimestamps1)
                spikeTimestamps1 = convertTimestamps( ts_raw, spikeTimestamps1, par.sr, 2 );
            end;
            
            %close(gcf);
            
            ts1 = [ts1 spikeTimestamps1];
            wavf1 = [wavf1;spikeWaveforms1];
            thr1 = [thr1 thr1];
            noise_std = [noise_std; noise_std_detect noise_std_sorted];
            
%             %% spike detection using osort
%             params_ts = [];
%             params_ts.nrNoiseTraces = 0;
%             params_ts.detectionMethod= 1;% 1 -> from power signal, 2 threshold positive, 3 threshold negative, 4 threshold abs, 5 wavelet
%             %params.detectionParams: depends on detectionMethod.
%             if params_ts.detectionMethod==1; params_ts.kernelSize = 18;end;
%             if params_ts.detectionMethod==4; params_ts.scaleRanges = [];params_ts.waveletName = '';end;
%             params_ts.alignMethod = 1;% 1-> find peak, 2->none, 3->peak of power signal, 4->peak of MTEO signal.
%             params_ts.extractionThreshold = 5;
%             params_ts.prewhiten = 0;
%             params_ts.samplingFreq = params.hdr.Fs;
%             params_ts.limit = 1000;
%             
%             [rawMean, filteredSignal, rawTraceSpikes,spikeWaveforms2, spikeTimestamps2, runStd2, upperlim, noiseTracesTmp] = extractSpikes(x', params.Hd, params_ts );
%             
%             wavf2o = spikeWaveforms2;%original waveform
%             
%             stdEstimates(zt) = std(filteredSignal);
%             
%             if size(noiseTracesTmp,1)>1
%                 noiseTracesTmp2=[];
%                 noiseTracesTmp2(1:size(noiseTracesTmp,1) ,1) = ones( size(noiseTracesTmp,1), 1 )*i;
%                 noiseTracesTmp2(1:size(noiseTracesTmp,1),2:size(noiseTracesTmp,2)+1)=noiseTracesTmp;
%                 noiseTraces = [noiseTraces; noiseTracesTmp2];
%             end
%             
%             %upsample and re-align
%             spikeWaveforms2=upsampleSpikes(spikeWaveforms2);
%             spikeWaveforms2 = realigneSpikes(spikeWaveforms2, spikeTimestamps2, params_ts.alignMethod, stdEstimates(zt));  %3==type is negative, do not remove any if too much shift
%             
%             %convert timestamps
%             if ~isempty(spikeTimestamps2)
%                 spikeTimestamps2 = convertTimestamps( ts_raw, spikeTimestamps2, params_ts.samplingFreq, 2 );
%             end;
%             
%             blockOffsets = [blockOffsets ts_raw(1)];
%             
%             wavf2 = [wavf2;spikeWaveforms2];
%             wavfo = [wavfo;wavf2o];
%             ts2 = [ts2 spikeTimestamps2];
%             thr2 = [ params_ts.extractionThreshold ];
            
            
        end;
        
        %% pool the spike indexes
        x = [];
        for it = 1:length(spkIdx)
            
            if it <2
                x = [x spkIdx{it}];
            else
                x = [x spkIdx{it}+sum([nsamp{it-1:-1:1}])];
            end;
        end;
        [spkIdx2] = x;
%         
%         %%
%         %whiten
%         if ( size(noiseTraces,1)>0 ) && ( size(wavfo,1)>0 )
%             [~, transUp, corr, ~] = posthocWhiten(noiseTraces, wavfo, params_ts.alignMethod);
%             wavfo = transUp; transUp = [];
%             
%             %only store the autocorrelation,not all noise traces
%             noiseTraces=corr; corr = [];
%         end;
%         
        %%        
        waveclus.spikes                     = wavf1;        
        waveclus.index                      = ts1;
        
        [dim] = size(wavf1);
        
        %%
        sortedSpikes1.newSpikeTimes     = [];
        sortedSpikes1.assignedCluster   = [];
        sortedSpikes1.wavf              = [];
        sortedSpikes1.num_clus          = [];
        
        if dim(1) > dim(2)% the number of spike events must larger than the number of samples in the waveform
            [sortedSpikes1] = doSpikeSorting_waveclus( waveclus , par );            
        end;
        waveclus = [];
        
%         %%
%         osort                               = [];
%         osort.ADbitVolts                    = params.hdr2.scalef;
%         osort.allSpikes                     = wavf2;    wavf2 = [];
%         osort.allSpikesCorrFree             = wavfo;
%         osort.allSpikesTimestamps           = ts2;
%         osort.blocksProcessed               = zt;
%         osort.noiseTraces                   = noiseTraces;
%         osort.blockOffsets                  = blockOffsets;
%         osort.dataSamplesRawUncorrected     = [];
%         osort.dataSamplesRaw                = dataSamples2;
%         osort.stdEstimates                  = stdEstimates;
%         osort.upperlim                      = upperlim;
%         osort.runStd2                       = runStd2;
%         osort.rawTraceSpikes                = rawTraceSpikes;
%         osort.filteredSignal                = filteredSignal;
%         osort.rawMean                       = rawMean;
%         
%         [sortedSpikes2] = doSpikeSorting_osort( osort );
        
        %%
        osort = [];
        rawMean = []; filteredSignal= []; rawTraceSpikes= [];spikeWaveforms2= []; spikeTimestamps2= []; runStd2= []; upperlim= []; noiseTracesTmp= [];
        
        %% make dummy ft structure for sorted spike data
        n = size(sortedSpikes1.wavf,2);
        wft = linspace(0,(n-1)/params.hdr.Fs,n);
        
        %timestamps defined with waveclus
        [spike1]                     =  struct;
        spike1.hdr                   = params.hdr;
        spike1.label                 = {['sig',sprintf('%03d',params.chanSel),'all_wvf']};
        spike1.timestamp{1}          = int64(sortedSpikes1.newSpikeTimes);
        spike1.waveform{1}(1,:,:)    = sortedSpikes1.wavf';%(:,1:2:end)';
        spike1.unit{1}               = sortedSpikes1.assignedCluster;
        spike1.params.hdr                   = params.hdr;
        spike1.dimord                = '{chan}_lead_time_spike';
        spike1.cfg.thr               = thr1;
        spike1.waveformtime          = wft;%(1:2:end);
        spike1.std                   = noise_std;
        sortedSpikes1 = [];
        
%         n = size(sortedSpikes2.newSpikesNegative,2);
%         wft = linspace(0,(n-1)/(params.hdr.Fs*4),n);
%         
%         %timestamps defined with O-sort
%         [spike2]                     =  struct;
%         spike2.hdr                   = params.hdr;
%         spike2.label                 = {['sig',sprintf('%03d',params.chanSel),'all_wvf']};
%         spike2.timestamp{1}          = int64(sortedSpikes2.newSpikesTimestampsNegative);
%         spike2.waveform{1}(1,:,:)    = sortedSpikes2.newSpikesNegative';%(:,1:8:end)';
%         spike2.unit{1}               = sortedSpikes2.assignedClusterNegative;
%         spike2.params.hdr                   = params.hdr;
%         spike2.dimord                = '{chan}_lead_time_spike';
%         spike2.cfg.thr               = thr2;
%         spike2.waveformtime          = wft;%(1:2:end);
%         spike2.std                   = stdEstimates;
%         
%         sortedSpikes2 = [];
        
        %% make dummy ft structure for unsorted spike data
        n = size(wavf1,2);
        wft = linspace(0,(n-1)/params.hdr.Fs,n);
        
        %timestamps defined with waveclus
        [spike3]                     =  struct;
        spike3.hdr                   = params.hdr;
        spike3.label                 = {['sig',sprintf('%03d',params.chanSel),'all_wvf']};
        spike3.timestamp{1}          = int64(ts1);
        spike3.waveform{1}(1,:,:)    = wavf1';%(:,1:2:end)';
        spike3.unit{1}               = zeros(1,length(ts1));
        spike3.params.hdr                   = params.hdr;
        spike3.dimord                = '{chan}_lead_time_spike';
        spike3.cfg.thr               = thr1;
        spike3.waveformtime          = wft;%(1:2:end);
        
        thr1 = [];
        ts1 = [];
        wavf1 = [];
        
%         n = size(wavfo,2);
%         wft = linspace(0,(n-1)/params.hdr.Fs,n);
%         
%         %timestamps defined with O-sort
%         [spike4]                     =  struct;
%         spike4.hdr                   = params.hdr;
%         spike4.label                 = {['sig',sprintf('%03d',params.chanSel),'all_wvf']};
%         spike4.timestamp{1}          = int64(ts2);
%         spike4.waveform{1}(1,:,:)    = wavfo';%(:,1:2:end)';
%         spike4.unit{1}               = zeros(1,length(ts2));
%         spike4.params.hdr                   = params.hdr;
%         spike4.dimord                = '{chan}_lead_time_spike';
%         spike4.cfg.thr               = thr2;
%         spike4.waveformtime          = wft;%(1:2:end);
%         
%         thr2 = [];
%         ts2 = [];
%         wavfo = [];
        
        %% do the spike time segmentation     
        [trl_dat] = make_spike_trials(params, spike1);
        cID = unique(trl_dat.unit);
        
        [spike_data1] = spike1;
        spike_data1.params                  = spike1.params;
        spike_data1.dimord                  = '{chan}_lead_time_spike';
        spike_data1.cfg                     = [];
        spike_data1.waveformtime            = spike1.waveformtime;
        spike_data1.trialtime               = trl_dat.trialtime;        
        for ct = 1:length(cID)
            sel_idx = find(trl_dat.unit == cID(ct));
            spike_data1.label{ct}                = ['sig',sprintf('%03d',ct),'sorted_wvf'];
            spike_data1.timestamp{ct}            = trl_dat.timestamp(sel_idx);
            spike_data1.waveform{ct}             = [];
            spike_data1.waveform{ct}(1,:,:)      = trl_dat.waveform(:,sel_idx);
            spike_data1.unit{ct}                 = trl_dat.unit(sel_idx);
            spike_data1.time{ct}                 = trl_dat.time(sel_idx);
            spike_data1.trial{ct}                = trl_dat.trial(sel_idx);
        end;
        
%         [trl_dat] = make_spike_trials(params, spike2);
%         cID = unique(trl_dat.unit);
%         
%         [spike_data2] = spike2;
%         spike_data2.params                  = spike2.params;
%         spike_data2.dimord                  = '{chan}_lead_time_spike';
%         spike_data2.cfg                     = [];
%         spike_data2.waveformtime            = spike2.waveformtime;
%         spike_data2.trialtime               = trl_dat.trialtime;
%         for ct = 1:length(cID)
%             sel_idx = find(trl_dat.unit == cID(ct));
%             spike_data2.time{ct}                 = trl_dat.time(sel_idx);
%             spike_data2.trial{ct}                = trl_dat.trial(sel_idx);
%             spike_data2.label{ct}                = ['sig',sprintf('%03d',ct),'sorted_wvf'];
%             spike_data2.timestamp{ct}            = trl_dat.timestamp(sel_idx);
%             spike_data2.waveform = [];
%             spike_data2.waveform{ct}(1,:,:)      = trl_dat.waveform(:,sel_idx);
%             spike_data2.unit{ct}                 = trl_dat.unit(sel_idx);
%         end;
        
        %cfg                         = [];
        %cfg.trl                     = [ double( params.TRL ) ];
        %cfg.trlunit                 = 'samples';
        %cfg.hdr                     = params.hdr;
        
        %[spike_data1]            = ft_spike_maketrials(cfg , spike1);
        %[spike_data2]            = ft_spike_maketrials(cfg , spike2);
        
        spike1 = [];
        spike2 = [];

        [trl_dat] = make_spike_trials(params, spike3);        

        [spike_data3] = spike3;
        spike_data3.label{1}                = ['sig',sprintf('%03d',params.chanSel),'all_wvf'];
        spike_data3.timestamp{1}            = trl_dat.timestamp;
        spike_data3.waveform                = [];
        spike_data3.waveform{1}(1,:,:)      = trl_dat.waveform;
        spike_data3.unit{1}                 = trl_dat.unit;
        spike_data3.params                  = spike3.params;
        spike_data3.dimord                  = '{chan}_lead_time_spike';
        spike_data3.cfg                     = [];        
        spike_data3.waveformtime            = spike3.waveformtime;
        spike_data3.time{1}                 = trl_dat.time;
        spike_data3.trial{1}                = trl_dat.trial;
        spike_data3.trialtime               = trl_dat.trialtime;
        
%         [trl_dat] = make_spike_trials(params, spike4);
%         
%         [spike_data4] = spike4;
%         spike_data4.label{1}                = ['sig',sprintf('%03d',params.chanSel),'all_wvf'];
%         spike_data4.timestamp{1}            = trl_dat.timestamp;
%         spike_data4.waveform                = [];
%         spike_data4.waveform{1}(1,:,:)      = trl_dat.waveform;
%         spike_data4.unit{1}                 = trl_dat.unit;
%         spike_data4.params                  = spike4.params;
%         spike_data4.dimord                  = '{chan}_lead_time_spike';
%         spike_data4.cfg                     = [];
%         spike_data4.waveformtime            = spike4.waveformtime;
%         spike_data4.time{1}                 = trl_dat.time;
%         spike_data4.trial{1}                = trl_dat.trial;
%         spike_data4.trialtime               = trl_dat.trialtime;
        
        %cfg                         = [];
        %cfg.trl                     = [ double( params.TRL ) ];
        %cfg.trlunit                 = 'samples';
        %cfg.hdr                     = params.hdr;
        
        %[spike_data3]            = ft_spike_maketrials(cfg , spike3);
        %[spike_data4]            = ft_spike_maketrials(cfg , spike4);
        spike3 = [];
        spike4 = [];
        
        %%
    end;
    
    if ( params.makeLFP ==1 ) || ( params.makeMUA ==1 )
        %% create dummy structure with CSC data
        
        % interpolate the timestamps to make time axis continuous
        [tsi] = timeStampinter(timestamps2);%interpolate timestamps
        tsi = tsi-min(tsi);
        tsi = tsi./1e6;
                
        [lfp_interp] = interpLFP(dataSamples,spkIdx2,[0.001 0.002],params.hdr2.Fs,'linear');

        %dummy ft structure
        CSC                         = [];
        CSC.hdr                     = params.hdr;
        CSC.fsample                 = params.hdr2.Fs;
        CSC.label                   = {chan};
        CSC.trial{1}                = lfp_interp;%dataSamples;%%
        CSC.time{1}                 = tsi;
        
        tsi = [];
        dataSamples = [];
                
        %% make the MUA
        if ( params.makeMUA ==1 )
            %bandpass spike-band
            cfg                     = [];
            cfg.bpfilter            = 'yes';
            cfg.bpfilttype          = 'fir';%'but';
            cfg.bpfreq              = [500 5000];
            %cfg.bpfiltord           = 2;
            cfg.padding             = 20;
            cfg.padtype             = 'mirror';
            
            [dum] = ft_preprocessing(cfg , CSC );
            
            %take the envelope
            cfg                     = [];
            cfg.hilbert             = 'abs';
            
            [dum] = ft_preprocessing(cfg, dum );
            
            % lowpass at 200 Hz
            cfg                     = [];
            cfg.lpfilter            = 'yes';
            cfg.lpfilttype          = 'fir';%'but';
            cfg.lpfreq              = 200;
            %cfg.lpfiltord           = 2;
            cfg.lpfilttype          = 'but';
            cfg.lpfiltord           = 2;
            
            [dum] = ft_preprocessing(cfg , dum );
            
            % segment trials
            cfg                     = [];
            cfg.trl                 = [ params.TRL ];
            
            [dum] = ft_redefinetrial(cfg, dum );
            
            %downsample
            cfg                     = [];
            cfg.resamplefs          = 1000;
            cfg.detrend             = 'no';
            
            [dum] = ft_resampledata(cfg , dum );
            
            %%            
            % notch filter line noise
            bsf = [ 50:50:200 ];
            for kt = 1:length(bsf)
                cfg                     = [];
                cfg.bsfilter            = 'yes';
                cfg.bsfreq              = [bsf(kt)-2 bsf(kt)+2];
                
                [dum] = ft_preprocessing(cfg , dum );
            end;
            
            [mua_data] = dum;
            dum = [];
            
        end;
        
        %% preprocessing of LFP data
        if (params.makeLFP ==1 )       
            
%             %%
%             cfg                     = [];
%             cfg.hpfilter            = 'yes';
%             cfg.hpfilttype          = 'but';%'fir';%
%             cfg.hpfiltord           = 4;
%             cfg.hpfreq              = 1;
%             %             cfg.padding             = 20;
%             %             cfg.padtype             = 'mirror';
%             
%             [CSC] = ft_preprocessing(cfg, CSC );
            
            cfg                     = [];
            cfg.lpfilter            = 'yes';
            cfg.lpfilttype          = 'but';%'fir';%
            cfg.lpfiltord           = 4;
            cfg.lpfreq              = 300;
            %             cfg.padding             = 20;
            %             cfg.padtype             = 'mirror';
            
            [CSC] = ft_preprocessing(cfg, CSC );
                        
            %%
            signal = CSC.trial;
            CSC.trial{1} = {};
            for zt = 1:length(CSC.label)
                for yt = 1:length(signal)
                    [cleanSignal,~] = CleanLineNoise( signal{yt}(zt,:) ,'Fs', CSC.fsample , 'noiseFreq', 50,'windowSize',100);
                    signal{yt}(zt,:) = cleanSignal;
                end;
                cleanSignal = [];
            end;
            CSC.trial = signal;
            signal = [];     
            
            %% do the LFP segmentation            
            try
                cfg                 = [];
                cfg.trl             = [ single(params.TRL) ];
                
                [CSC] = ft_redefinetrial(cfg, CSC );
            catch
                cfg                 = [];
                cfg.trl             = [ params.TRL ];
                
                [CSC] = ft_redefinetrial(cfg, CSC );
            end;                   

            %%
            cfg                         = [];
            cfg.detrend                 = 'no';
            cfg.resamplefs              = 1000;
            
            [CSC]                  = ft_resampledata(cfg , CSC );
            
            %%
            [lfp_data] = CSC;
            CSC = [];
        end;
    end;
    
catch ME
    rethrow(ME)
end;

%%
if ( params.makeSpikes == 1 )
    data{end+1} = spike_data1;%sorted WC
    %data{end+1} = spike_data2;%sorted OS
    data{end+1} = spike_data3;%unsorted WC
    %data{end+1} = spike_data4;%unsorted OS
    
    readme{end+1} = 'sortedSpikeTimesWC';
    %readme{end+1} = 'sortedSpikeTimesOsort';
    readme{end+1} = 'unsortedSpikeTimesWC';
    %readme{end+1} = 'unsortedSpikeTimesOsort';
end;

if ( params.makeLFP == 1 )
    data{end+1} = lfp_data;
    readme{end+1} = 'LFPdata';
end;

if ( params.makeMUA == 1 )
    data{end+1} = mua_data;
    readme{end+1} = 'MUAdata';
end;

return;
