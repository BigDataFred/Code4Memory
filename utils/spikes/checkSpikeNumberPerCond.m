function [chck] = checkSpikeNumberPerCond(c1Idx,c2Idx,tIx,data_spk,Ntrl,mode)
% input: spike count matrix (values must be binary) amd indexes of
% condition trials

if isempty(c1Idx) && isempty(c2Idx) && isempty(tIx)
    
    chck = zeros(1,length(data_spk));
    for it = 1:length(data_spk)
        x = data_spk{it};
        cnt = sum(x(:));% count total number of spikes for condition1
        if cnt > Ntrl
            chck(it) = 1;
        end;
    end;

elseif isempty(c1Idx) && isempty(c2Idx) && ~isempty(tIx)
    
    chck = zeros(1,length(data_spk));
    for it = 1:length(data_spk)
        x = data_spk{it}(tIx,:);
        cnt = sum(x(:));% count total number of spikes for condition1
        if cnt > Ntrl
            chck(it) = 1;
        end;
    end;
    
else
    
    chck = zeros(1,length(data_spk));
    
    for it = 1:length(data_spk)% loop over the channels
        
        x1 = data_spk{it}(:,c1Idx);
        x2 = data_spk{it}(:,c2Idx);
        
        cnt1 = sum(x1(:));% count total number of spikes for condition1
        cnt2 = sum(x2(:));% count total number of spikes for condition2
        
        if mode == 1
            if cnt1 > Ntrl && cnt2 > Ntrl
                chck(it) = 1;% passed the minimum of spikes in both conditions
            end;
        else
            if (cnt1 + cnt2) > Ntrl
                chck(it) = 1;% passed the minimum of spikes in both conditions
            end;
        end;
    end;
    
end;

return
