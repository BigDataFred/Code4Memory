function [spike_dat] = segment_spike_data(et,st,params)
%%
  
[pSTH] = compute_pSTH(et,st,[params.pre params.post],params.bw);
[trlst] = st2trl(st,et,params.Fs,[params.pre params.post]);

[spike_dat] = struct;
spike_dat.dt = linspace(-1,1,size(pSTH,2));

spike_dat.psth = pSTH;
spike_dat.trlst = trlst;
spike_dat.frate = compute_firing_rate(spike_dat.psth,params.bw);
spike_dat.frate = smoothFR(spike_dat.frate,params.bw*params.Fs,params.window);
spike_dat.units = 's';

return;