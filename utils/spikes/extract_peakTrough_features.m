function [dt,pt] = extract_peakTrough_features(wavf,time_vec)

Y = squeeze(wavf);

dt = [];
pt = [];
for it = 1:size(Y,2)
    
    [v1,ix1] = max(Y(:,it));
    [v2,ix2] = min(Y(:,it));
    
    dt(it) = diff(time_vec(sort([ix1 ix2],'ascend')));
    pt(it) = v2-v1;
    
end;