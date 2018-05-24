function [spkTms,cond,dt] = clusterSpikeTimes2mat(selIx1,selIx2,spikeDat,trlIdx,ix,toi,bw)
%% extract spike times to matrix
spkTms = cell(1,length(selIx1));
cond = cell(1,length(selIx1));

for it = 1:length( selIx1)
    
    % get the trials of the spike times
    trl = spikeDat{selIx1(it)}.trial{selIx2(it)+1};
    trl = trl(find(ismember(trl,trlIdx)));% retain trials during learning
            
    % get each trials ID
    trlId = trlIdx;
    n = length(trlId);% total number of trials in cluster
    
    % get the spike times
    ts = spikeDat{selIx1(it)}.time{selIx2(it)+1}.*1e3;%convert to ms
    
    x =cell(1,n);    
    % trials with cond1 - cond2
    for jt = 1:n;
        ix2 = find(trl == trlId(jt));   
        if isempty(ix2)
            disp('toto');
        end;
        x{jt} = ts(ix2);
    end;    
    
    if ~isempty(ix)
        for ft = 1:length(ix)
            if size(ix{ft},1) > size(ix{ft},2)
                ix{ft} = ix{ft}';
            end,
        end;
        
        [cond{it}] = [ix{1} ix{2} ix{3};...
            ones(1,length(ix{1})) 2*ones(1,length(ix{2})) 3*ones(1,length(ix{3}))];
    end;
    
    dt = toi(1):bw:toi(2);
    spkTms{it} = zeros(n,length(dt));
    for jt = 1:length(x);
        
        x2 = x{jt};
        x2(x2<dt(1)) = [];
        x2(x2>dt(end)) = [];
        
        spkTms{it}(jt,:) =hist(x2,dt);
        spkTms{it}(jt,:) = conv(spkTms{it}(jt,:),gausswin(10),'same');
        
    end;
    
end;