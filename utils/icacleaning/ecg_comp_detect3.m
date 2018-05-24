function [ecg_indx] = ecg_comp_detect3(comp)

%% detection of ICs with ECG artifacts
% transform IC time-courses into zscores

ts= zeros(length(comp.time),1);for it = 1:length(comp.time);ts(it) = length(comp.time{it});end;ts = sum(ts); 
%%
concat = zeros(length(comp.label),ts);
%%
indx = 1:length(comp.time{1});
for it = 1:length(comp.trial)
    
    
    concat(:,indx) =comp.trial{it};
    if it <length(comp.trial)
        indx = indx(end)+1:indx(end)+length(comp.time{it+1});
    end;
end;
%%
ecg_indx = [];
k = 0;
zsc = (concat - repmat(mean(concat,2),[1 size(concat,2)]))./repmat(std(concat,1,2),[1 size(concat,2)]);
zsc = abs(zsc) >= 4-k;

while k<2
    zsc = zsc(setdiff(1:size(zsc,1),ecg_indx),:);
    
    
    for jt = 1:size(zsc,1)
        
        cross_indx = [];
        sample_distance = [];
        del_indx = [];
        
        cross_indx = find(zsc(jt,:));
        sample_distance = diff(cross_indx);
        del_indx = find(sample_distance<comp.fsample*0.5);
        zsc(jt,cross_indx(del_indx)) = 0;
        
    end;
    
    minutes = size(zsc,2)/(comp.fsample*60);
    ecg_freq_min = 40*minutes;% minimum frequency of 45Hz
    ecg_freq_max = 130*minutes;% maximum frequency of 150Hz
   
    % keep index of IC that matches physiological heartbeat frequency
    [ecg_indx] = [ecg_indx;find(sum(zsc,2) >= ecg_freq_min & sum(zsc,2) <= ecg_freq_max)];
    
    if ~isempty(ecg_indx)
        k = k+1;
    else
        error('the data has no ECG-component');
    end;
end;
%return