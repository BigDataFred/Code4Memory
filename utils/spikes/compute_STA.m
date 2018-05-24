function [STA] = compute_STA(win,n,lfp)            

STA = zeros(win*2+1,size(n,2));
for jt = 1:size(n,2)%loop over trials
    
    spk_ix = find(n(:,jt)~=0);% find spike time indexes
    
    for kt = 1:length(spk_ix)%loop over number of spikes
        if (spk_ix(kt)-win >0) && (spk_ix(kt)+win <=size(n,1))
            x = lfp(spk_ix(kt)-win:spk_ix(kt)+win,jt);
            STA(:,jt) = STA(:,jt) + x./length(spk_ix);%divide by n-spikes
        end;
    end;
    
end;
