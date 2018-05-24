function [emg_idx] = emg_comp_detect(comp_freq)

%%
d = zeros(length(comp_freq.label),1);
for it = 1:length(comp_freq.label)
    
    Y = squeeze(mean(comp_freq.powspctrm(:,it,:),1));
    Y = (Y-min(Y))./(max(Y)-min(Y));
    
    X = comp_freq.freq';
    X =X(find(sign(comp_freq.freq-30)==1 & sign(comp_freq.freq-max(comp_freq.freq))==-1));
    Y= Y(find(sign(comp_freq.freq-30)==1 & sign(comp_freq.freq-max(comp_freq.freq))==-1));
    f = comp_freq.freq(find(sign(comp_freq.freq-max(comp_freq.freq))==1 & sign(comp_freq.freq-max(comp_freq.freq))==-1));
    
    Y(find(sign(f-47)==1 & sign(f-53)==-1))= NaN;
    
    xi = find(isnan(Y));
    Y(xi) = [];
    X(xi) = [];
    
    [b,bint,r,rint,stats] = regress(Y,[ones(size(X)) X]);
    
    yp = b(1)+b(2).*X;
    
    d(it) = sum((yp-Y).^2);
    
end;
d = d-mean(d);
d = d./std(d);
[d,s_idx]  = sort(d);
%%
[emg_idx] = sort(s_idx(find(sign(d-0.1) ==1)));
% %%
% %K = zeros(1,length(comp_freq.label));%kurtosis
% %S = zeros(1,length(comp_freq.label));%skewness
% P = zeros(1,length(comp_freq.label));%power
% for it = 1:length(comp_freq.label)
%     
% %     K(it) = median(kurtosis(squeeze(comp_freq.powspctrm(:,it,:)),0,2));
% %     S(it) = median(skewness(squeeze(comp_freq.powspctrm(:,it,:)),0,2));
%     
%     x = squeeze((comp_freq.powspctrm(:,it,:)));
% %    x = (x-min(min(x)))./(max(max(x))-min(min(x)));
%     
%     %x = (x-repmat(mean(x,2),[1 size(x,2)]));
%     %x = x./repmat(std(x,0,2),[1 size(x,2)]);
%     x = cumsum(x,2);
%     x = max(x(:,comp_freq.freq>=30),[],2);
%     P(it) = median(x); 
% end;
% %%
% % K = (K - mean(K))./std(K);
% % S =(S - mean(S))./std(S);
% %%
% trshP = (P)-(mean(P));
% trshP = trshP./std(trshP);
% trshP = (trshP)>2;
% %%
% [v,six] = sort(sum(trshP,1));
% %v  = diff(v)./diff(1:length(v));
% %v = diff(v)./diff(1:length(v));
% v = (v-mean(v))./std(v);
% %%
% E = zeros(size(comp_freq.powspctrm,1),length(comp_freq.label));%entropy
% for jt = 1:size(comp_freq.powspctrm,1)
%     for it = 1:length(comp_freq.label)
%         
%         x1 = squeeze(comp_freq.powspctrm(jt,it,:));
%         x1 = x1./sum(x1);
%         %x1 = -sum(log(x1).*x1);
%         
%         x2 = squeeze(mean(comp_freq.powspctrm(jt,setdiff(1:length(comp_freq.label),find(trshP==1)),:),2))';
%         x2 = x2./sum(x2);
%         %x2 = -sum(log(x2).*x2);
%         
%         E(jt,it) = sum(x2.*log(x2./x1'));%(x1-x2)./x2;
%     end;
% end;
% 
% e = median(E,1);
% e = (e-mean(e))./std(e);
%%
% % idx1 = find(K<1);
% % idx2 = find(S<1.5);
% idx3 = six(min(find(v>2)):length(v));
% idx4 = find(e<1);
% 
% %idx5 = intersect(idx1,idx2);
% idx6 = intersect(idx3,idx4);
% %comp_idx5 = intersect(idx5,idx6);
% comp_idx5 = idx6;
% emg_idx = comp_idx5;

% %% visualize
%     %%
%     figure;
%     set(gcf,'Position',[50 130 1730 820]);
%     subplot(2,length(emg_idx)+3,1);
%     imagesc(P(:,six));
%     title('Sorted 20-80Hz power [\sigma]')
%     xlabel('Component #');
%     ylabel('Sinle trial #');
% 
%     subplot(2,length(emg_idx)+3,2);
%     [ax,h1,h2] = plotyy(1:length(v),v,1:length(e),e(six));
%     set(h1,'LineStyle','.');
%     set(h2,'LineStyle','.');
% 
%     hold(ax(1),'on');
%     hold(ax(2),'on');
%     plot(ax(1),[1 length(v)],[2 2],'b--');
%     plot(ax(2),[1 length(e)],[1 1],'g--');
%     axis(ax,'tight');
%     xlabel(ax(1),'Component #');
%     ylabel(ax(1),'Median 20-80Hz power [\sigma]');
%     ylabel(ax(2),'Median KL-divergence [\sigma]');
% 
%     subplot(2,length(emg_idx)+3,3);
%     [ax,h1,h2] = plotyy(1:length(k),k,1:length(s),s(six));
%     set(h1,'LineStyle','.');
%     set(h2,'LineStyle','.');
%     hold(ax(1),'on');
%     hold(ax(2),'on');
%     plot(ax(1),[1 length(k)],[5 5],'b--');
%     plot(ax(2),[1 length(s)],[2 2],'g--');
%     axis(ax,'tight');
% 
%     xlabel(ax(1),'Component #');
%     ylabel(ax(1),'Median kurtosis');
%     ylabel(ax(2),'Median skewness');
% 
%     subplot(2,length(emg_idx)+3,(length(emg_idx)+3)+1);
%     hold on;
%     plot(k,s,'.','Color',[.75 .75 .75]);
%     if ~isempty(comp_idx5)
%         plot(k(comp_idx5),s(comp_idx5),'r.');
%     end;
%     plot([5 5],[min(s) max(s)],'r--');
%     plot([min(k) max(k)],[2 2],'r--');
%     xlabel('Kurtosis');
%     ylabel('Skewness');
%     axis tight;
% 
%     subplot(2,length(emg_idx)+3,(length(emg_idx)+3)+2);
%     hold on;
%     plot(e,v,'.','Color',[.75 .75 .75]);
%     v = sum(trshP,1);
%     v = (v-mean(v))./std(v);
% 
%     if ~isempty(comp_idx5)
%         plot(e(comp_idx5),v(comp_idx5),'r.');
%     end;
%     plot([1 1],[min(v) max(v)],'r--');
%     plot([min(e) max(e)],[2 2],'r--');
%     xlabel('KL-Divergence [\sigma]');
%     ylabel('20-80Hz Power [\sigma]');
%     axis tight;