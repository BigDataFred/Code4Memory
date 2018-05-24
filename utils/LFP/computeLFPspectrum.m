function [spectrumLFP] = computeLFPspectrum(data_lfp,Fs)

TW = 5;
params                  = [];
params.tapers           = [TW , 2*TW-1];
params.pad              = 2;
params.Fs               = Fs;
params.err              = 0;%[1 0.001];
params.trialave         = 0;
params.fpass            = [1 30];


%[S1,fx1,SErr1] = mtspectrumc(data_lfp,params);
[S1,fx1] = mtspectrumc(data_lfp,params);

% TW = floor(single(size(data_lfp,1)/Fs)*10);
% params                  = [];
% params.tapers           = [TW , 2*TW-1];
% params.pad              = 2;
% params.Fs               = Fs;
% params.err              = 0;%[1 0.001];
% params.trialave         = 1;
% params.fpass            = [30 200];
% 
% %[S2,fx2,SErr2] = mtspectrumc(data_lfp,params);
% [S2,fx2] = mtspectrumc(data_lfp,params);

spectrumLFP = [];
spectrumLFP.fx1 = fx1;
% spectrumLFP.fx2 = fx2;
spectrumLFP.S1 = S1;
% spectrumLFP.S2 = S2;
% spectrumLFP.SErr1 = SErr1;
% spectrumLFP.SErr2 = SErr2;