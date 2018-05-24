function [data,readme,chan] = analyze_spontaneous_data(p2NLX,NLXsesh,CSCidx,makeSpikes,makeLFP,makeMUA)

%%
for it = 1:length(CSCidx)
    
    sehsPath = [p2NLX,filesep,NLXsesh,filesep];
    
    CSC_files = dir([sehsPath,'CSC_*.ncs']);
    
    %% read the CSC data
    FieldSelection(1) = 1;%timestamps
    FieldSelection(2) = 0;
    FieldSelection(3) = 0;%sample freq
    FieldSelection(4) = 0;
    FieldSelection(5) = 1;%samples
    ExtractHeader = 1;
    
    ExtractMode = 1; % 2 = extract record index range; 4 = extract timestamps range.
    
    [timestamps, dataSamples,hdr] = Nlx2MatCSC_v3([sehsPath,CSC_files(CSCidx(it)).name], FieldSelection, ExtractHeader, ExtractMode, []);
    
    chck = regexp(hdr,'ADBitVolts');
    selIdx = [];
    for jt = 1:length(chck);selIdx(jt) = ~isempty(chck{jt});end;
    selIdx = find(selIdx~=0);
    scalef = str2double(hdr{selIdx}(min(regexp(hdr{selIdx},'\d')):end));
    
    
    chck = regexp(hdr,'SamplingFrequency');
    selIdx = [];
    for jt = 1:length(chck);selIdx(jt) = ~isempty(chck{jt});end;
    selIdx = find(selIdx~=0);
    Fs = str2double(hdr{selIdx}(min(regexp(hdr{selIdx},'\d')):end));
    
    %% flatten
    dataSamples=double(dataSamples(:))';
    dataSamples = dataSamples.*scalef.*1e6;
    
     %% estimate threshold parameters for clustering
     [par]                               = set_parameters_Bham( Fs );
     
     [~,chan,~] = fileparts( [sehsPath,CSC_files(CSCidx(it)).name] );
     
     par.fnamespc                        = [ par.fnamespc,chan ] ;
     par.fname_in                        = [ par.fname_in,chan ];
     
     [fs,~,~,~] = amp_detect(dataSamples,par);
     
     waveclus                            = [];
     waveclus.thr                        = par.sort_fmin * median(abs(fs))/0.6745;
     waveclus.thrmax                     = par.stdmax * median(abs(fs))/0.6745;
    
    %% make the spike times
    nsegs = 6;
    
    lts = length(timestamps);
    segl = floor(lts/nsegs);
    
    tsmin = 1: segl : lts;
    tsmin = tsmin(1: nsegs );
    tsmax = tsmin - 1;
    tsmax = tsmax(2:end);
    tsmax = [tsmax, lts];
    
    tsmin = timestamps(int64(tsmin));
    tsmax = timestamps(int64(tsmax));
    
     %% make the spike data
    if makeSpikes ==1
        
        %% pre-allocate
        ts1 = [];        thr1 = [];        wavf1 = []; noise_std =[];
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
            
            [ ts_raw, dataSamples2 ] = Nlx2MatCSC_v3([sehsPath,CSC_files(CSCidx(it)).name], FieldSelection, ExtractHeader, ExtractMode, [fromInd toInd]);
         
            x = double(dataSamples2(:))';
            x = x*scalef*1e6;
            
            %% spike detection using wave_clus
             [~,spikeWaveforms1,thr1,spikeTimestamps1,noise_std_detect,noise_std_sorted] = amp_detect(x,par);
            %spikeTimestamps1 = spikeTimestamps1.*1e6./params.hdr.Fs+tsmin(zt);%convert samples to timestamps and add start time of system clock
            
             if ~isempty(spikeTimestamps1)
                spikeTimestamps1 = convertTimestamps( ts_raw, spikeTimestamps1, par.sr, 2 );
            end;
            
            close(gcf);
            
            ts1 = [ts1 spikeTimestamps1];
            wavf1 = [wavf1;spikeWaveforms1];
            thr1 = [thr1 thr1];
            noise_std = [noise_std; noise_std_detect noise_std_sorted];                                  
            
        end;                
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
        
        %% make dummy ft structure for unsorted spike data
        n = size(sortedSpikes1.wavf,2);
        wft = linspace(0,(n-1)/Fs,n);
        
        %timestamps defined with waveclus
        [spike1]                     =  struct;
        spike1.hdr                   = hdr;
        spike1.label                 = {['sig_',chan,'_all_wvf']};
        spike1.timestamp{1}          = int64(sortedSpikes1.newSpikeTimes);
        spike1.waveform{1}(1,:,:)    = sortedSpikes1.wavf';%(:,1:2:end)';
        spike1.unit{1}               = sortedSpikes1.assignedCluster;
        spike1.params.hdr                   = hdr;
        spike1.dimord                = '{chan}_lead_time_spike';
        spike1.cfg.thr               = thr1;
        spike1.waveformtime          = wft;%(1:2:end);
        spike1.std                   = noise_std;
        sortedSpikes1 = [];
        
        n = size(wavf1,2);
        wft = linspace(0,(n-1)/Fs,n);
        
        %timestamps defined with waveclus
        [spike2]                     =  struct;
        spike2.hdr                   = hdr;
        spike2.label                 = {['sig_',chan,'_all_wvf']};
        spike2.timestamp{1}          = int64(ts1);
        spike2.waveform{1}(1,:,:)    = wavf1';%(:,1:2:end)';
        spike2.unit{1}               = zeros(1,length(ts1));
        spike2.params.hdr                   = hdr;
        spike2.dimord                = '{chan}_lead_time_spike';
        spike2.cfg.thr               = thr1;
        spike2.waveformtime          = wft;%(1:2:end);
        
        thr1 = [];
        ts1 = [];
        wavf1 = [];
        
        [spike_data1] = spike1;
        spike_data1.dimord                  = '{chan}_lead_time_spike';
        spike_data1.cfg                     = [];
        spike_data1.waveformtime            = spike1.waveformtime;
        if~isempty(spike_data1.unit{:})
            cID = unique(spike_data1.unit{:});
            for ct = 1:length(cID)
                sel_idx = find(spike1.unit{:} == cID(ct));
                spike_data1.label{ct}                = ['sig',sprintf('%03d',ct),'sorted_wvf'];
                spike_data1.timestamp{ct}            = spike1.timestamp{1}(sel_idx);
                spike_data1.waveform{ct}             = [];
                spike_data1.waveform{ct}(1,:,:)      = spike1.waveform{1}(:,:,sel_idx);
                spike_data1.unit{ct}                 = spike1.unit{1}(sel_idx);
            end;
        end;
        
        [spike_data2] = spike2;
        spike_data2.dimord                  = '{chan}_lead_time_spike';
        spike_data2.cfg                     = [];
        spike_data2.waveformtime            = spike2.waveformtime;
        if~isempty(spike_data2.unit{:})
            cID = unique(spike_data2.unit{:});
            for ct = 1:length(cID)
                sel_idx = find(spike2.unit{:} == cID(ct));
                spike_data2.label{ct}                = ['sig',sprintf('%03d',ct),'sorted_wvf'];
                spike_data2.timestamp{ct}            = spike2.timestamp{1}(sel_idx);
                spike_data2.waveform{ct}             = [];
                spike_data2.waveform{ct}(1,:,:)      = spike2.waveform{1}(:,:,sel_idx);
                spike_data2.unit{ct}                 = spike2.unit{1}(sel_idx);
            end;
        end;
        %%
    end;

        if ( makeLFP ==1 ) || ( makeMUA ==1 )
        %% create dummy structure with CSC data
        
        % interpolate the timestamps to make time axis continuous
        [tsi] = timeStampinter(timestamps);%interpolate timestamps
        tsi = tsi-min(tsi);
        tsi = tsi./1e6;
        
        %dummy ft structure
        CSC                         = [];
        CSC.hdr                     = hdr;
        CSC.fsample                 = Fs;
        CSC.label                   = {chan};
        CSC.trial{1}                = dataSamples;
        CSC.time{1}                 = tsi;
        
        tsi = [];
        dataSamples = [];
                
        %% make the MUA
        if (makeMUA ==1 )
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
        if (makeLFP ==1 )                        
            
            %%
            signal = CSC.trial;
            CSC.trial{1} = {};
            for zt = 1:length(CSC.label)
                for yt = 1:length(signal)
                    [cleanSignal,~] = CleanLineNoise( signal{yt}(zt,:) ,'Fs', CSC.fsample , 'noiseFreq', 50,'windowSize',1);
                    signal{yt}(zt,:) = cleanSignal;
                end;
                cleanSignal = [];
            end;
            CSC.trial = signal;
            signal = [];
            
            %%
            cfg                     = [];      
            cfg.lpfilter            = 'yes';
            cfg.lpfilttype          = 'but';
            %cfg.lpfiltord           = 2;
            cfg.padding             = 20;
            cfg.padtype             = 'mirror';            
            cfg.lpfreq              =  300;
            
            [CSC] = ft_preprocessing(cfg, CSC );

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
        
end;

%%
data = {};
readme = {};
if ( makeSpikes == 1 )
    data{end+1} = spike_data1;%sorted WC
    data{end+1} = spike_data2;%unsorted WC
    
    readme{end+1} = 'sortedSpikeTimesWC';
    readme{end+1} = 'unsortedSpikeTimesWC';
end;

if ( makeLFP == 1 )
    data{end+1} = lfp_data;
    readme{end+1} = 'LFPdata';
end;

if ( makeMUA == 1 )
    data{end+1} = mua_data;
    readme{end+1} = 'MUAdata';
end;
