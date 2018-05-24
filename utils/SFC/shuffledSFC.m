function [pval,Cemp,Cshuf] = shuffledSFC(data_lfp,data_spk,params,Fs)

if nargin <3   
    Fs = 1e3;
    TW = (size(data_lfp,1)/Fs)*1/(size(data_lfp,1)/Fs);
    params                  = [];
    params.tapers           = [TW , 2*TW-1];
    params.pad              = -1;
    params.Fs               = Fs;
    params.err              = 0;
    params.trialave         = 1;
    params.fpass            = [0 30];
end;
    
if params.pad ~=-1
    fx = params.fpass(1):1/(2^(nextpow2(size(data_lfp,1))+(params.pad))/Fs):params.fpass(2);
else
    fx = params.fpass(1):1/(size(data_lfp,1)/Fs):params.fpass(2);
end;

[Cemp]=coherencycpb(data_lfp,data_spk,params);

Cshuf = zeros(1e3,length(fx));
parfor it = 1:1e4
    fprintf([num2str(it),'/',num2str(1e3)]);
    trl = 1:size(data_lfp,2);
    trl = trl(randperm(length(trl)));
    
    dumLFP = data_lfp;
    dumLFP = dumLFP(:,trl);
    
    [Cshuf(it,:)]=coherencycpb(dumLFP,data_spk,params);    
    fprintf('\n');
end;

pval = zeros(size(Cshuf,2),1);
parfor it = 1:size(Cshuf,2)
    pval(it) = (length(find(Cshuf(:,it)>=Cemp(it)))/size(Cshuf,1));        
end;