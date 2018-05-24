function [sortedSpikes] = doSpikeSorting_osort( osortdata )
%%
filesToProcess = 1;

paramsIn = [];
paramsIn.rawFileVersion = 2; %1 is analog cheetah, 2 is digital cheetah (NCS), 3 is txt file.  determines sampling freq&dynamic range.
paramsIn.samplingFreq = osortdata; %only used if rawFileVersion==3

%which tasks to execute
paramsIn.tillBlocks = 999;  %how many blocks to process (each ~20s). 0=no limit.
paramsIn.doDetection = 0;
paramsIn.doSorting = 1;
paramsIn.doFigures = 0;
paramsIn.noProjectionTest = 1;
paramsIn.doRawGraphs = 0;
paramsIn.doGroundNormalization=0;

paramsIn.displayFigures = 1 ;  %1 yes (keep open), 0 no (export and close immediately); for production use, use 0 
paramsIn.minNrSpikes=50; %min nr spikes assigned to a cluster for it to be valid
                                                                                                                         
%params
paramsIn.blockNrRawFig=[ 1 2 ];
%paramsIn.outputFormat='png';
paramsIn.thresholdMethod=1; %1=approx, 2=exact
paramsIn.prewhiten=0; %0=no, 1=yes,whiten raw signal (dont)
paramsIn.defaultAlignMethod=3;  %only used if peak finding method is "findPeak". 1=max, 2=min, 3=mixed
paramsIn.peakAlignMethod=1; %1 find Peak, 2 none, 3 peak of power, 4 MTEO peak
                        
%for power detection method
paramsIn.detectionMethod=1; %1 power, 2 T pos, 3 T min, 3 T abs, 4 wavelet
dp.kernelSize=18; 
paramsIn.detectionParams=dp;
extractionThreshold = 5;  % extraction threshold

thres         = [repmat(extractionThreshold, 1, length(filesToProcess))];
%%
handles = [];

[samplingFreq, limit, postfix] = defineFileFormat(paramsIn.rawFileVersion, paramsIn.samplingFreq)
handles.samplingFreq = samplingFreq; %sampling freq of raw data
handles.limit = limit; %dynamic range

handles.blocksProcessedForInit                      = osortdata.blocksProcessed;
handles.dataSamplesRaw                              = osortdata.dataSamplesRaw;
handles.rawMean                                     = osortdata.rawMean;
handles.rawTraceSpikes                              = osortdata.rawTraceSpikes;
handles.runStd2                                     = osortdata.runStd2;
handles.upperlim                                    = osortdata.upperlim;
handles.filteredSignal                              = osortdata.filteredSignal;
handles.noiseTraces                                 = osortdata.noiseTraces;

%estimate s.d. of raw signal
handles.stdEstimateOrig = calculateStdEstimate(osortdata.stdEstimates); %mean(stdEstimates);

handles.allSpikesNegative                           = osortdata.allSpikes;
handles.allSpikesTimestampsNegative                 = osortdata.allSpikesTimestamps;
handles.newSpikesNegative                           = osortdata.allSpikes;
handles.newSpikesTimestampsNegative                 = osortdata.allSpikesTimestamps;
handles.spikesSolvedNegative                        = osortdata.allSpikes;
handles.allSpikesNoiseFree                          = [];
handles.allSpikesCorrFree                           = osortdata.allSpikesCorrFree;

%for compatibility reasons
handles.allSpikesPositive                           =[];
handles.allSpikesTimestampsPositive                 =[];
handles.newSpikesPositive                           =[];
handles.newSpikesTimestampsPositive                 =[];
handles.spikesSolvedPositive                        =[];

%New: also store the scaling factor
handles.scalingFactor                               = osortdata.ADbitVolts;
handles.correctionFactorThreshold                   =0; 
handles.stdEstimate                                 = handles.stdEstimateOrig*handles.correctionFactorThreshold;

handles.paramExtractionThreshold                    = thres;

handles = copyStructFields( paramsIn, handles, {'minNrSpikes','blockNrRawFig','doGroundNormalization','rawFileVersion','detectionMethod','detectionParams','peakAlignMethod','displayFigures'});

handles.normalizationChannels                       = [];
handles.doGroundNormalization                       = 0;
handles.normalizationChannels                       = [];

[handles] = sortMain( [], handles, 2, paramsIn.thresholdMethod  ); %2=no GUI

sortedSpikes = handles;