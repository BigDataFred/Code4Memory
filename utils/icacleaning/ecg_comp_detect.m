function [ecg_indx] = ecg_comp_detect(comp)
%% detects ICs with ECG artifacts
%
%
%% initialize
ecg_indx = [];
%% concatenate ICA data
concat = zeros(length(comp.label),length(comp.time{1})*length(comp.trial));
indx = 1:length(comp.time{1});
for it = 1:length(comp.trial)        
    concat(:,indx) =comp.trial{it};
    if it <length(comp.trial)
        indx = indx(end)+1:indx(end)+length(comp.time{it+1});
    end;
end;
%z-score IC-data
zsc = (concat - repmat(mean(concat,2),[1 size(concat,2)]))./repmat(std(concat,1,2),[1 size(concat,2)]);
clear concat;
%%
k = 0;% counter
for it = 1:4
    
    zsc2 = [];
    zsc2 = abs(zsc) >= it;
    
    for jt = 1:size(zsc2,1)
        
        cross_indx = [];
        sample_distance = [];
        del_indx = [];
        
        cross_indx = find(zsc2(jt,:));
        sample_distance = diff(cross_indx);
        del_indx = find(sample_distance<comp.fsample*0.25);
        zsc2(jt,cross_indx(del_indx)) = 0;
        
    end;
    
    minutes = size(zsc2,2)/(comp.fsample*60);
    ecg_freq_min = 25*minutes;% minimum frequency of 45Hz
    ecg_freq_max = 190*minutes;% maximum frequency of 150Hz
    
    dum = sum(zsc2,2);
    
    [v,idx] = sort(dum);
    v = (v-mean(v))./std(v);

    dum = dum.*(dum>=ecg_freq_min);
    dum = dum.*(dum<=ecg_freq_max);
    
    % keep index of IC that matches physiological heartbeat frequency    
    ecg_indx = [ecg_indx;intersect(idx(find(v>2)),find(dum~=0))];
    %ecg_indx = [ecg_indx;idx(find(v>2))];
end;
ecg_indx = unique(ecg_indx);
return