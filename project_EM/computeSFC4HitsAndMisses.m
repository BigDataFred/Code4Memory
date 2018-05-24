function [staDat] = computeSFC4HitsAndMisses(data_lfp,data_lfpBP,data_spk,hitIdx,missIdx,timeGrid,Fs)

%% SFC = STA/STP
rand('state',sum(100*clock));

smp =timeGrid;
D = [-.5 .5];
plt = 'n';
err = 0;
T = [timeGrid(1) timeGrid(end)];

[allSTA]    = sta(data_spk,data_lfp,smp,plt,[],T,D,err);
[allSTAbp]  = sta(data_spk,data_lfpBP,smp,plt,[],T,D,err);

[powSTA]    = computeSpectrumAndSFC4timeRange({mean(allSTA,1)'}, [], [], [],Fs);
[allSTP]    = computeSpectrumAndSFC4timeRange({allSTA'}, [], [], [],Fs);

[allSFC] = (powSTA{1}.S1./mean(allSTP{1}.S1,2)).*100;

[hitsSTA]   = sta(data_spk(hitIdx), data_lfp(:,hitIdx), smp,plt,[],T,D,err);
[missSTA]   = sta(data_spk(missIdx), data_lfp(:,missIdx), smp,plt,[],T,D,err);
[hitsSTAbp] = sta(data_spk(hitIdx), data_lfpBP(:,hitIdx), smp,plt,[],T,D,err);
[missSTAbp] = sta(data_spk(missIdx), data_lfpBP(:,missIdx),smp,plt,[],T,D,err);

[n] = [size(hitsSTA,1) size(missSTA,1)];
[~,ix] = min(n);
% if ix ~=2
%     error('Misses have more spikes than Hits');
% end;
if ix == 2
    
    [powSTA] = computeSpectrumAndSFC4timeRange({mean(missSTA,1)'}, [], [], [],Fs);
    [missSTP] = computeSpectrumAndSFC4timeRange({missSTA'}, [], [], [],Fs);
        
    missSFC = (powSTA{1}.S1./mean(missSTP{1}.S1,2)).*100;%MISSES
    
    nIter = 500;
    x = zeros(nIter,length(missSFC));
    y = zeros(nIter,length(missSTP{1}.S1));
    
    parfor zt = 1:nIter

        randSel = randperm(max(n));
        randSel = randSel(1:min(n));
        
        dum1 = {hitsSTA(randSel,:)'};
        
        [dumSTA1] = computeSpectrumAndSFC4timeRange( {mean(dum1{1},2)}, [], [], [],Fs);
        [dumSTP1] = computeSpectrumAndSFC4timeRange(dum1, [], [], [],Fs);
        
        x(zt,:) = (dumSTA1{1}.S1./mean(dumSTP1{1}.S1,2)).*100;
        y(zt,:) = dumSTP1{1}.S1;
        
    end;
    hitsSFC = mean(x,1);%HITS
    hitsSTP = mean(y,1);
else
    [powSTA] = computeSpectrumAndSFC4timeRange({mean(hitsSTA,1)'}, [], [], [],Fs);
    [hitsSTP] = computeSpectrumAndSFC4timeRange({hitsSTA'}, [], [], [],Fs);
    
    
    hitsSFC = (powSTA{1}.S1./mean(hitsSTP{1}.S1,2)).*100;%MISSES
    
    nIter = 500;
    x = zeros(nIter,length(hitsSFC));
    y = zeros(nIter,length(hitsSTP{1}.S1));
    
    parfor zt = 1:nIter
        
        randSel = randperm(max(n));
        randSel = randSel(1:min(n));
        
        dum1 = {missSTA(randSel,:)'};
        
        [dumSTA1] = computeSpectrumAndSFC4timeRange( {mean(dum1{1},2)}, [], [], [],Fs);
        [dumSTP1] = computeSpectrumAndSFC4timeRange(dum1, [], [], [],Fs);
        
        x(zt,:) = (dumSTA1{1}.S1./mean(dumSTP1{1}.S1,2)).*100;
        y(zt,:) = dumSTP1{1}.S1;
        
    end;
    missSFC = mean(x,1);%HITS
    missSTP = mean(y,1);
end;

%%
staDat.allSTA      = allSTA;
staDat.allSTAbp    = allSTAbp;
staDat.allSTP      = allSTP;
staDat.allSFC      = allSFC;
staDat.hitsSFC     = hitsSFC;
staDat.missSFC     = missSFC;
staDat.hitsSTA     = hitsSTA;
staDat.missSTA     = missSTA;
staDat.hitsSTP     = hitsSTP;
staDat.missSTP     = missSTP;
staDat.hitsSTAbp   = hitsSTAbp;
staDat.missSTAbp   = missSTAbp;

