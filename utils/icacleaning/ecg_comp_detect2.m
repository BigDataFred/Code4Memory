function [ecg_indx] = ecg_comp_detect2(cfg,comp_data)
%% detects ICs with ECG artifacts
%
%
%% define physiological range of heartbeat frequency
if isfield(cfg,'ecg_freq_min');else cfg.ecg_freq_min = 30;end;% = 0.5 heartbeat per 1 sec
if isfield(cfg,'ecg_freq_max');else cfg.ecg_freq_max = 180;end;% = 3 heartbeats per 1 sec
if isfield(cfg,'treshold');else cfg.treshold = 2;end;% treshold in sd
%% concatenate segmented ICA data into continuous time series
concat = zeros(length(comp_data.label),length(comp_data.time{1})*length(comp_data.trial));
indx = 1:length(comp_data.time{1});
for it = 1:length(comp_data.trial)        
    concat(:,indx) =comp_data.trial{it};
    if it <length(comp_data.trial)
        indx = indx(end)+1:indx(end)+length(comp_data.time{it+1});
    end;
end;
%% z-score IC-data
zsc = (concat - repmat(mean(concat,2),[1 size(concat,2)]))./repmat(std(concat,1,2),[1 size(concat,2)]);
clear concat;
if isfield(cfg,'ecg_indx')
    zsc = zsc(cfg.ecg_indx,:);
end;
%% transform to absolute values to facilitate treshold detection
%%
zsc1 = zsc.*(zsc >= cfg.treshold);%do first-pass tresholding
zsc2 = zsc.*(zsc <= -cfg.treshold);%do first-pass tresholding
f = zeros(1,size(zsc,1));
for jt = 1:size(zsc,1)
    
    if sum(zsc1(jt,:) >= cfg.treshold,2) > sum(zsc2(jt,:)<= -cfg.treshold)
        x = zsc1(jt,:);
    else
        x = abs(zsc2(jt,:));
    end;
    tx = 0:1/comp_data.fsample:(size(x,2)-1)/comp_data.fsample;
    
    x = x > median(x(x~=0))/1.5;%do median-based second-pass tresholding
    idx = find(sign(diff(x))==1)+1;%find zero-crossings
    dx = diff(idx);% inter-heart-beat distance in samples
    
    del_idx = [find(sign(dx./comp_data.fsample-.5)==-1) find(sign(dx./comp_data.fsample-1.5)==1)];
    
    dx(del_idx) = [];
    idx(del_idx) = [];
    if (length(idx)/(tx(end)/60))/60 >.65
        if length(find(dx >=median(dx)-round(median(dx)*0.125) & dx <= median(dx)+round(median(dx)*0.125)))/length(dx) >.3% check the minimum number of occurrences of inter-heart-beat intervals
            if length(idx) >= cfg.ecg_freq_min*(tx(end)/60) && length(idx) <= cfg.ecg_freq_max*(tx(end)/60)
                f(jt) = 1;
            end;
        end;
    else
        f(jt) = 0;
    end;
end;
%% keep index of IC that matches physiological heartbeat frequency
% and get the indexes of ICs carrying EKG activity
%[ecg_indx] = find(f<cfg.ecg_freq_max & f>cfg.ecg_freq_min)';
[ecg_indx] = find(f==1)';
%return