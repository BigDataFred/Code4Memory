function [paramsTF1,paramsTF2,movingwin1,movingwin2,paramsSP1,paramsSP2] = setParameters4spectralAnalysis(Fs,tw,trialave)

%% parameters for time-frequency analysis (I)

%lower frequencies
movingwin1 = [.8 .2];
T = movingwin1(1);% length of time window in s
W = 2;% smoothing (+/- 1Hz)
TW = T*W;% time-bandwidth product
k = floor(2*(T*W)-1);% number of tapers

paramsTF1              = [];
paramsTF1.Fs           = Fs;
paramsTF1.pad          = 2;
paramsTF1.fpass        = [0.5 30];
paramsTF1.tapers       = [TW k];
if trialave == false
    paramsTF1.trialave     = 1;
else
    paramsTF1.trialave     = 0;
end;
%higher frequencies
movingwin2 = [0.25 0.0625];
T = movingwin2(1);% length of time window in s
W = 8;% smoothing (+/- 10Hz)cd
TW = T*W;% time-bandwidth product
k = floor(2*(T*W)-1);% number of tapers

paramsTF2              = [];
paramsTF2.Fs           = Fs;
paramsTF2.pad          = 4;
paramsTF2.fpass        = [20 220];
paramsTF2.tapers       = [TW k];
if trialave == false
    paramsTF2.trialave     = 1;
else
    paramsTF2.trialave     = 0;
end;

%% parameters for spectral analysis (II)
T = tw;
W = 1;%smoothing (+/- 1Hz)
TW = T*W;
k = 2*TW-1;

paramsSP1              = [];
paramsSP1.Fs           = Fs;
paramsSP1.pad          = 2;
paramsSP1.fpass        = [0.5 30];
paramsSP1.tapers       = [TW k];
paramsSP1.trialave     = 1;
paramsSP1.err          = [2 0.001];%jacknife error

T = tw;
W = 4;%smoothing (+/- 4Hz)
TW = T*W;
k = 2*TW-1;

paramsSP2              = [];
paramsSP2.Fs           = Fs;
paramsSP2.pad          = 2;
paramsSP2.fpass        = [20 220];
paramsSP2.tapers       = [TW k];
paramsSP2.trialave     = 1;
paramsSP2.err          = [2 0.001];%jacknife error
