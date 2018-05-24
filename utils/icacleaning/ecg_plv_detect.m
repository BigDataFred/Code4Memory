function [ecg_idx] = ecg_plv_detect(cfg,comp_freq)
%this section of code uses the ICs detected during the first-pass
%as "seed"-components to look for additional ICs carrying EKG
%activity

if ~isfield(cfg,'ecg_idx');cfg.ecg_idx = [];end;

if ~isempty(cfg.ecg_idx)
    %% measure phase locking between seed EKG components and all other ICs
    cfg2 = [];
    cfg2.method = 'plv';
    cfg2.channelcmb = [];
    add_ekg_label = cell(length(cfg.ecg_idx),1);
    for it = 1:length(cfg.ecg_idx)
                
        cfg2.channelcmb = [repmat(comp_freq.label(cfg.ecg_idx(it)),[length(comp_freq.label)-length(cfg.ecg_idx) 1]) comp_freq.label(setdiff(1:length(comp_freq.label),cfg.ecg_idx))];
        alt_idx = find(~ismember(1:length(cfg.ecg_idx),it));                
        
        [ecg_coh] = ft_connectivityanalysis(cfg2,comp_freq);
        %% z-score the PLV values across ICs
        for kt = 1:size(ecg_coh.plvspctrm,2)
            ecg_coh.plvspctrm(:,kt) = (ecg_coh.plvspctrm(:,kt)-mean(ecg_coh.plvspctrm(:,kt)))./(std(ecg_coh.plvspctrm(:,kt)));
        end;
        
        dum = ecg_coh.plvspctrm>10;%treshold the z-scored PLV values
        dum = sum(dum,2);%calculate the sum of the surface below the PLV spectrum
        dum = (dum-mean(dum))./std(dum);%z-score of the sum under the curve
        
        sel_idx = unique(sort([find(dum>4)]));% find components which have large sum
        add_ekg_label{it} = unique(sort(ecg_coh.labelcmb(sel_idx,2)));%get the label informatino of the selected ICs
    end
    %use the label information to get to the index referring to the
    %ICs which show a strong degree of coherence with the EKG
    %component detected during the first-pass
    k = 0;
    ecg_idx = zeros(length(add_ekg_label),1);
    for kt = 1:length(add_ekg_label)
        for nt = 1:length(add_ekg_label{kt})
        if ~isempty(find(strcmp(add_ekg_label{kt}(nt),comp_freq.label(:))))
            k = k+1;
            ecg_idx(k) = find(strcmp(add_ekg_label{kt}(nt),comp_freq.label(:)));
        end;
        end;
    end;
    ecg_idx(k+1:end) = [];
    
    if size(ecg_idx,1)>size(ecg_idx,2)
        ecg_idx = ecg_idx';
        ecg_idx = unique(sort(ecg_idx));
    end;
else
    ecg_idx = [];
end;