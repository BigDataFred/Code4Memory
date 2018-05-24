function [spkDat] = extractClusterDat(spkDat,uIx)

spkDat = spkDat{uIx(1)};

spkDat.label                    = spkDat.label(uIx(2));
spkDat.timestamp                = spkDat.timestamp(uIx(2));
spkDat.waveform                 = spkDat.waveform(uIx(2));
spkDat.unit                     = spkDat.unit(uIx(2));
spkDat.time                     = spkDat.time(uIx(2));
spkDat.trial                    = spkDat.trial(uIx(2));