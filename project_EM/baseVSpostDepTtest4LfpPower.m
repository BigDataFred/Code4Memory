function [statInfo] = baseVSpostDepTtest4LfpPower(lfpDat,foi,toi,chanIx)

%% recruit workers for parallel computing
if isempty(gcp('nocreate'))
    parpool(36,'SpmdEnabled', false)
end;
    
%%
if isempty(chanIx)
    [nMWs] = length( lfpDat.LFPseg );
    chanIx = 1:nMWs;
end;

%%
d = toi(2)-toi(1);

T = d;
W = 3;

paramsLF                  =[];
paramsLF.Fs               = lfpDat.Fs;
paramsLF.fpass            = [foi(1) foi(2)];
paramsLF.pad              = 5;
paramsLF.tapers           = [T*W  2*T*W-1];
paramsLF.trialave         = 0;

[nMWs] = length( chanIx );

powBase = cell(nMWs,1);
powPost = cell(nMWs,1);
empT = zeros(nMWs,1);
parfor curMW = 1:length( chanIx )
    fprintf([num2str(curMW),'/',num2str(nMWs)]);
    
    [ntrl] = size(lfpDat.LFPseg{chanIx(curMW)},2);
    if ntrl >=20
        [S1,~] = mtspectrumc( gradient( lfpDat.LFPseg{chanIx(curMW)}(lfpDat.trlTime>=-d & lfpDat.trlTime <=0,:)' )', paramsLF);
        [S2,~] = mtspectrumc( gradient( lfpDat.LFPseg{chanIx(curMW)}(lfpDat.trlTime>=toi(1) & lfpDat.trlTime <=toi(2),:)' )', paramsLF);
        
        powBase{curMW} = mean(20*log10(S1),1);
        powPost{curMW} = mean(20*log10(S2),1);
        rtn= sqrt(length(powPost{curMW}));
        empT(curMW) = mean(powPost{curMW}- powBase{curMW})./std(powPost{curMW}- powBase{curMW})*rtn;
    else
        powBase{curMW} = [];
        powPost{curMW} = [];        
    end;
    
    fprintf('\n');
end;

%%
nPerm = 2000;
permT = zeros(nMWs,nPerm);
for curMW = 1:nMWs
    if ~isempty(powBase{curMW}) && ~isempty(powPost{curMW})
        parfor curPer = 1:nPerm
            fprintf([num2str(curPer),'/',num2str(nPerm)]);
            x = [powPost{curMW}' powBase{curMW}'];
            rIx = randperm(size(x,1));
            rIx = rIx(1:ceil(length(rIx)/2));
            tmp = x(rIx,:);
            x(rIx,1) = tmp(:,2);
            x(rIx,2) = tmp(:,1);
            rtn= sqrt(size(x,1));
            permT(curMW,curPer) = mean(x(:,2)- x(:,1))./std(x(:,2)- x(:,1))*rtn;
            fprintf('\n');
        end;
    end;
end;

%%
pval = zeros(nMWs,2);
for it = 1:size(permT,1)
    [pval(it,1)] = sum(permT(it,:)>= empT(it))/nPerm;
    [pval(it,2)] = sum(permT(it,:)<= empT(it))/nPerm;
end;

[selIxP] = chanIx(find(pval(:,1)<2.5e-2));
[selIxN] = chanIx(find(pval(:,2)<2.5e-2));

if size(selIxP,2)>size(selIxP,1)
    selIxP = selIxP';    
end;

if size(selIxN,2)>size(selIxN,1)
    selIxN = selIxN';    
end;

%%
statInfo.pval = pval;
statInfo.selIxP = selIxP;
statInfo.selIxN = selIxN;
statInfo.permT = permT;
statInfo.empT = empT;


