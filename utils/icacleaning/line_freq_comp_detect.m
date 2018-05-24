function [line_idx] = line_freq_comp_detect(comp_freq)
%%

pf = zeros(size(comp_freq.powspctrm,1),length(comp_freq.label));
for jt = 1:size(comp_freq.powspctrm,1) 
    for it = 1:length(comp_freq.label)
       
        x = squeeze(mean(comp_freq.powspctrm(jt,it,:),1));
        pf(jt,it) = comp_freq.freq(find(x==max(x)));
    end;
end;

pf1 = (pf >= 16 & pf <=18);
pf2 = (pf >= 48 & pf <=52);

Spf = sum(pf1,1)+sum(pf2,1);
Spf = Spf/size(comp_freq.powspctrm,1);

[line_idx] = unique(find(Spf >=.75));
%%

% sel_idx = median(find(comp_freq.freq >= 49.9 & comp_freq.freq <=50.1));
% figure;
% for it = 1:length(line_idx)
%     subplot(6,2,it)
%     plot(squeeze(comp_freq.powspctrm(:,line_idx(it),sel_idx)),'b.');
%     axis tight;
%     title(Spf(line_idx(it)));
% end;
% 
% sel_idx = median(find(comp_freq.freq >= 16 & comp_freq.freq <=18));
% figure;
% for it = 1:length(line_idx)
%     subplot(6,2,it)
%     plot(squeeze(comp_freq.powspctrm(:,line_idx(it),sel_idx)),'b.');
%     axis tight;
%     title(Spf(line_idx(it)));
% end;