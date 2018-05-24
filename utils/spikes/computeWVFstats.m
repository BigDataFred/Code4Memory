function [wvfStats] = computeWVFstats(x)

%%
wvfStats = [];
[n,m] = size(x);
if n == 64
    tmp = mean(x,2);
    r = zeros(size(x,2),1);
    parfor it = 1:size(x,2)
        r(it) = corr(x(:,it),tmp,'Type','Spearman');
    end;
    
    [pix] = find(sign(r)==1);
    
    ppct = length(pix)/length(r);
    m=tmp;
    sd = std(x,0,2);
    
    [pa,~] = max(m);
    [ta,~] = min(m);
    p2p = diff([pa ta]);
    
    res = [];
    for it = 1:size(x,2)
        res = [res;x(:,it)-m];
    end;
    snr = abs(p2p)/(std(res)*5);
    
    [~,ix1] = max(abs(m));
    dx = gradient(m);
    if ix1+1 == length(dx)
        ix2 = ix1;
    else
        try
            if (sign(m(ix1))==1) && (sign(dx(ix1+1)) ==-1)
                ix2 = ix1+min(find(sign(dx(ix1+1:end))==1));
            elseif (sign(m(ix1))==-1) && (sign(dx(ix1+1)) ==1)
                ix2 = ix1+min(find(sign(dx(ix1+1:end))==-1));
                
            else
                error('cannot determine waveform-duration');
            end;
        catch
            return;
        end;
    end;
    
    [~,Fy] = gradient(x);% Fx = derivative over trials (horizontal), Fy = derivative over time (vertical)
    
    wvfStats.m = m;
    wvfStats.sd = sd;
    wvfStats.snr =snr;
    wvfStats.r = r;
    wvfStats.rpct = ppct;
    wvfStats.peakAmp = pa;
    wvfStats.peakIx = ix1;
    wvfStats.troughAmp = ta;
    wvfStats.troughIx = ix2;
    wvfStats.wvf = x;
    wvfStats.peakDur = (round(abs(diff([mean(ix1) mean(ix2)])))/32000).*1e3;
    wvfStats.peak2trough = diff([pa ta],[],2);
end;
return;