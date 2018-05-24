function [emg_idx2] = emg_coh_detect(emg_idx,comp_freq)

emg_idx2 = [];

cfg = [];
cfg.method = 'coh';
cfg.channelcmb = [];
for it = 1:length(emg_idx)
    
    cfg.channelcmb = [repmat(comp_freq.label(emg_idx(it)),[1 length(comp_freq.label)-length(emg_idx)])' comp_freq.label(setdiff(1:length(comp_freq.label),emg_idx))];
    
    emg_coh = ft_connectivityanalysis(cfg,comp_freq);
    if strcmp(cfg.method,'coh')
        for jt = 1:size(emg_coh.cohspctrm,2)
            emg_coh.cohspctrm(:,jt) = emg_coh.cohspctrm(:,jt) - repmat(mean(emg_coh.cohspctrm(:,jt),1),[length(emg_coh.labelcmb) 1]);
            emg_coh.cohspctrm(:,jt) = emg_coh.cohspctrm(:,jt) ./ repmat(std(emg_coh.cohspctrm(:,jt),0,1),[length(emg_coh.labelcmb) 1]);
        end;
        emg_coh.cohspctrm = emg_coh.cohspctrm.*(emg_coh.cohspctrm>2);
        
        y = squeeze(sum(emg_coh.cohspctrm,2));
    else
        for jt = 1:size(emg_coh.plvspctrm,2)
            emg_coh.plvspctrm(:,jt) = emg_coh.plvspctrm(:,jt) - repmat(mean(emg_coh.plvspctrm(:,jt),1),[length(emg_coh.labelcmb) 1]);
            emg_coh.plvspctrm(:,jt) = emg_coh.plvspctrm(:,jt) - repmat(std(emg_coh.plvspctrm(:,jt),0,1),[length(emg_coh.labelcmb) 1]);
        end;
        y = squeeze(sum(emg_coh.plvspctrm,2));
    end;
    y = y-mean(y);
    y = y./std(y);
    y = y.*(y>2);
    [y,s_idx] = sort(y);
    ydx = diff(y);
    ydx = ydx-mean(ydx);
    ydx = ydx./std(ydx);
    
    emg_idx2 = sort(unique([emg_idx2;s_idx([min(find(ydx>=4)):length(ydx)]-1)]));
end;