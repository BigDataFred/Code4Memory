function [rms_idx] = return_rms_idx(comp)
%%
cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = [30 100];

[bf_comp] = ft_preprocessing(cfg,comp);

rms = zeros(length(bf_comp.trial),length(bf_comp.label));
for it = 1:length(bf_comp.trial)
    
        rms(it,:) = mean(sqrt(abs(bf_comp.trial{it}).^2),2);
        
end;

rms = std(rms,0,1);

rms = rms-mean(rms);
rms = rms./std(rms);

[rms_idx] = find(sign(rms-2)==1);