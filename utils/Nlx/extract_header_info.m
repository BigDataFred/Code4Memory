function [hdr] = extract_header_info( filename )
 
%%
[timestamps,nrBlocks,nrSamples,sampleFreq,isContinous,headerInfo] = getRawCSCTimestamps( filename );
%%
x = zeros(length(headerInfo),1);
for it = 1:length(headerInfo)
    x(it) = ~isempty(regexp(headerInfo{it},'-ADBitVolts'));
end;
scale_factor = str2double(headerInfo{find(x==1)}(regexp(headerInfo{find(x==1)},'-ADBitVolts')+11:end));
%%
x = zeros(length(headerInfo),1);
for it = 1:length(headerInfo)
    x(it) = ~isempty(regexp(headerInfo{it},'-AcqEntName'));
end;
chan_label = headerInfo{find(x==1)}(regexp(headerInfo{find(x==1)},'-AcqEntName')+12:end);
%%
hdr.nChans = length(chan_label);
hdr.label = {chan_label};
hdr.filenmae = filename;
hdr.nTrials = [];
hdr.Fs = sampleFreq;
hdr.nSamplesPre = [];
hdr.nSamples = nrSamples;
hdr.FirstTimeStamp = min(sort(timestamps));
hdr.LastTimeStamp = max(sort(timestamps));
hdr.TimeStampPerSample = (1/sampleFreq)*1e6;
hdr.orig = [];
hdr.chantype = {};
hdr.chanunit = {};
hdr.ADBitVolts = scale_factor;

  