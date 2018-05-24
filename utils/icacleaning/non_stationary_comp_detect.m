
function [nst_idx] = non_stationary_comp_detect(comp)
%% select the variability across trials of components
SD = zeros(length(comp.label),length(comp.trial));
for it = 1:length(comp.label)
    for jt = 1:length(comp.trial)
        SD(it,jt) = var(abs(hilbert(comp.trial{jt}(it,:))).^2);
    end;
    SD(it,:) = SD(it,:) - mean(SD(it,:),2);
    SD(it,:) = SD(it,:)./std(SD(it,:),0,2);
end;
SD2 = abs(SD)>=2;
%%
sd = sum(SD2,2)./size(SD2,2);
sd = squeeze(std((SD),0,2));
sd = (sd-mean(sd))./std(sd,0,1);
[sd,s_idx] = sort(sd);
roc = ((diff(sd)))./diff(1:length(sd))';
roc = (roc-mean(roc))./(std(roc));

idx = min(intersect(find(sd > 2),find(roc >2)+1));
idx = idx:length(sd);
nst_idx = unique(sort(s_idx(idx)))';

