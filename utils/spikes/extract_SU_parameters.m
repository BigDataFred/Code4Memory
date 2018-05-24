function [SU] = extract_SU_parameters(spk_dat,sel_ix)


SU.wvf = squeeze( spk_dat.waveform{sel_ix} )';
SU.ts = spk_dat.time{sel_ix}.*1e3;% time in ms

% uses a routine from Ueli Ruthishaeuser
[SU2.wvf,SU2.ts, ~] = postDetectionFilter( SU.wvf, SU.ts);

SU.sd = mean(spk_dat.std(:,2));
SU.t = spk_dat.waveformtime;

return;