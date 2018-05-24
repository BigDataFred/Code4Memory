function [empT1,empT2,pval1,pval2,thr1,thr2] = extractSigRespsFR(FR,nperms)

%% compute STATS
empT1 = [];    empT2 = []; % stats
for ct = 1:length(FR)
    [~,~,~,stat1] = ttest(FR{ct}(:,2),FR{ct}(:,1));% compare spike count cue vs base
    empT1(ct) = stat1.tstat;% save empirical t-stat
    
    [~,~,~,stat2] = ttest(FR{ct}(:,3),FR{ct}(:,1));% compare spike count encoding vs base
    empT2(ct) = stat2.tstat;% save empirical t-stat
end;

%% create null-distribution from permutation data
randT1 = [];        randT2 = [];
randF1 = [];        randF2 = [];

parfor pt = 1:nperms
    %for pt = 1:nperms
    
    fprintf([num2str(pt),'/',num2str(nperms)]);
    
    rT1 = [];        rT2 = [];
    rF1 = [];        rF2 = [];
    
    for it = 1:length(FR)%loop over all clusters
        
        x = [];    x = FR{it}(:,[1 2]);% base and cue
        x2 = [];
        for kt = 1:size(x,1)%loop over trials
            x2(kt,:) = x(kt,randperm(size(x,2)));  % randomly swap the condition label
        end;
        
        [~,~,~,stat1] = ttest(x2(:,2),x2(:,1));% compare spike count cue vs base
        rT1(it) = stat1.tstat;% save empirical t-stat
        
        x =[];    x = FR{it}(:,[1 3]);% base and encoding
        x2 = [];
        
        for kt = 1:size(x,1)%loop over trials
            x2(kt,:) = x(kt,randperm(size(x,2))); % randomly swap the condition label
        end;
        
        [~,~,~,stat2] = ttest(x2(:,2),x2(:,1));% compare spike count encoding vs base
        rT2(it) = stat2.tstat;% save empirical t-stat
        
        
    end;% end of loop across clusters
    
    % keep the stats from this specific random partition
    randT1(pt,:) = [min(rT1) max(rT1)];
    randT2(pt,:) = [min(rT2) max(rT2)];
    
    fprintf('\n');
    
end;% end of loop across permutations

%% compute p-values
pval1 = zeros(length(empT1),1);
for it = 1:length(empT1);
    if isnan(empT1(it));
        pval1(it) = 1;
    else
        %if sign(empT1(it))==1
            pval1(it) = length(find(randT1(:,2)>=empT1(it)))/nperms;
        %else
        %    pval1(it) = length(find(randT1(:,1)<=empT1(it)))/nperms;
        %end;
    end;
end;

pval2 = zeros(length(empT2),1);
for it = 1:length(empT2);
    if isnan(empT2(it));
        pval2(it) = 1;
    else
        %if sign(empT2(it))==1
            pval2(it) = length(find(randT2(:,2)>=empT2(it)))/nperms;
        %else
        %    pval2(it) = length(find(randT2(:,1)<=empT2(it)))/nperms;
        %end;
    end;
end;


%%
thr1 = zeros(length(FR),1);
thr1(find(pval1 <= 0.025))  = 1;
%thr1(intersect(find(pval1 < 0.025),find(sign(empT1)==1)))  = 1;
%thr1(intersect(find(pval1 < 0.025),find(sign(empT1)==-1))) = -1;

thr2 = zeros(length(FR),1);
thr2(find(pval2 <= 0.025))  = 1;
%thr2(intersect(find(pval2 < 0.025),find(sign(empT2)==1)))  = 1;
%thr2(intersect(find(pval2 < 0.025),find(sign(empT2)==-1))) = -1;
