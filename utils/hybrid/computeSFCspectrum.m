function [spectrumLFP] = computeSFCspectrum(data_lfp, data_spk, Fs)

TW = 5;
params                  = [];
params.tapers           = [TW , 2*TW-1];
params.pad              = 2;
params.Fs               = Fs;
params.err              = 0;%[1 0.001];
params.trialave         = 1;
params.fpass            = [1 30];


%[SFC1,~,~,~,~,fx1,~,~,~,CErr1] = coherencycpb(data_lfp, data_spk, params);
[SFC1,~,~,~,~,fx1] = coherencycpb(data_lfp, data_spk, params);

% TW = floor(single(size(data_lfp,1)/Fs)*10);
% params                  = [];
% params.tapers           = [TW , 2*TW-1];
% params.pad              = 2;
% params.Fs               = Fs;
% params.err              = 0;%[1 0.001];
% params.trialave         = 1;
% params.fpass            = [30 200];
% 
% %[SFC2,~,~,~,~,fx2,~,~,~,CErr2] = coherencycpb(data_lfp, data_spk, params);
% [SFC2,~,~,~,~,fx2] = coherencycpb(data_lfp, data_spk, params);

spectrumLFP = [];
spectrumLFP.fx1 = fx1;
% spectrumLFP.fx2 = fx2;
spectrumLFP.SFC1 = SFC1;
% spectrumLFP.SFC2 = SFC2;
% spectrumLFP.SFC1 = phi1;
% spectrumLFP.SFC2 = phi2;
% spectrumLFP.SFC1 = confC1;
% spectrumLFP.SFC2 = confC2;
% spectrumLFP.SFC1 = phistd1;
% spectrumLFP.SFC2 = phistd2;
% spectrumLFP.CErr1 = CErr1;
% spectrumLFP.CErr2 = CErr2;