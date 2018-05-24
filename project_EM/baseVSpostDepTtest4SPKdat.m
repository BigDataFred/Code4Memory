function [statInfo] = baseVSpostDepTtest4SPKdat(baseCnt,postCnt)

%% recruit workers for parallel computing
if isempty(gcp('nocreate'))
    parpool(36,'SpmdEnabled', false)
end;

%%
rtn= sqrt(size(postCnt,1));
empT = mean(postCnt - baseCnt)./std(postCnt-baseCnt)*rtn;

%%
nPerm = 2000;
permT = zeros(nPerm,1);
if ~isempty(baseCnt) && ~isempty(postCnt)
    parfor curPer = 1:nPerm
        fprintf([num2str(curPer),'/',num2str(nPerm)]);
        x = [baseCnt postCnt];
        rIx = randperm(size(x,1));
        rIx = rIx(1:ceil(length(rIx)/2));
        tmp = x(rIx,:);
        x(rIx,1) = tmp(:,2);
        x(rIx,2) = tmp(:,1);        
        permT(curPer) = mean(x(:,2)- x(:,1))./std(x(:,2)- x(:,1))*rtn;
        fprintf('\n');
    end;
end;

%%
[pval(1)] = sum(permT>= empT)/nPerm;
[pval(2)] = sum(permT<= empT)/nPerm;


%%
statInfo.pval = pval;
statInfo.permT = permT;
statInfo.empT = empT;

