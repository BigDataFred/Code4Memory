function [ecg_idx,qrsA] = clean_ecg(ica_data)
%%
addpath('~froux/froux/sl_meg/mcode/icacleaning/');
%% do a first pass to detect EKG-components
cfg =[];
cfg.ecg_freq_min = 30;
cfg.ecg_freq_max = 210;
cfg.treshold = 2;

[ecg_idx1] = ecg_comp_detect2(cfg,ica_data);
%% do second pass to detect wether first pass returned genuine EKG-components
cfg = [];
cfg.ecg_idx = ecg_idx1;

[ecg_idx1,qrsA] = ft_component_class_ekg2(cfg,ica_data);
%%
if size(ecg_idx1,2)>size(ecg_idx1,1)
    ecg_idx1 = ecg_idx1';
end;
%% do third pass to detect components that are phase-locked to ECG-signal
cfg = [];
cfg.method = 'mtmfft';
cfg.taper = 'dpss';
cfg.tapsmofrq = round(1/max(ica_data.time{1})*1000)/1000*10/2;
cfg.pad = 'maxperlen';
cfg.foi = 0:200;
cfg.keeptrials = 'yes';
cfg.output = 'fourier';

[comp_cf] = ft_freqanalysis(cfg,ica_data);

cfg = [];
cfg.ecg_idx = ecg_idx1;

[ecg_idx2] = ecg_plv_detect(cfg,comp_cf);
%%
if size(ecg_idx2,2)>size(ecg_idx2,1)
    ecg_idx2 = ecg_idx2';
end;
%%
if size(ecg_idx1,1) ~= size(ecg_idx2,1)
    
    if size(ecg_idx1,1) > size(ecg_idx1,2)
        ecg_idx1 = ecg_idx1';
    end;

    if size(ecg_idx2,1) > size(ecg_idx2,2)
        ecg_idx2 = ecg_idx2';
    end;

end;
[ecg_idx] = sort([ecg_idx1,ecg_idx2]);