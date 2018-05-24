function [comp_data_idx] = eog_comp_detect2(comp_data,eeg_data)
%% low pass 

cfg = [];
cfg.lpfilter = 'yes';
cfg.lpfreq = 15;

[eeg_data] = ft_preprocessing(cfg,eeg_data);
[comp_data] = ft_preprocessing(cfg,comp_data);
%% search for correlation between EOG-signals and ICs in the time domain

r1 = zeros(length(comp_data.trial),1);
r2 = zeros(length(comp_data.trial),1);
x1 = zeros(length(comp_data.label),1);
x2 = zeros(length(comp_data.label),1);
for it = 1:length(comp_data.label)
    for jt = 1:length(comp_data.trial)
        [r1(jt)] = max(abs(xcorr(comp_data.trial{jt}(it,:)',eeg_data.trial{jt}(1,:)','coeff')));
        [r2(jt)] = max(abs(xcorr(comp_data.trial{jt}(it,:)',eeg_data.trial{jt}(2,:)','coeff')));
    end;
    % correlation between EEG and IC time series needs to be high
    x1(it) =  length(find(abs(r1)>.45))/length(r1);
    x2(it) =  length(find(abs(r2)>.45))/length(r2);
end;
%% artefact needs to be present in more than 50% of trials
eog_idx = [find(x1>.5);find(x2>.5)];
%% return index of ICs that show high similarity with EOG-signals
comp_data_idx = sort(unique(eog_idx));