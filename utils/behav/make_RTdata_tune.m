function [RTdat] = make_RTdata_tune(LogDat)
RTdat.cond = LogDat.dat(2:end,3);
RTdat.stimID = str2double(LogDat.dat(2:end,1));
RTdat.resp = str2double(LogDat.dat(2:end,4));
RTdat.RT = str2double(LogDat.dat(2:end,end));
%%
M = nanmean(RTdat.RT);
SD = nanstd(RTdat.RT);
z = (RTdat.RT-M)./SD;
RTdat.oL = find(z >= 3.5);

[v,s_idx] = sort(RTdat.RT);

RTdat.s_idx = s_idx;