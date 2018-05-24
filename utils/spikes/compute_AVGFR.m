function [FR] = compute_AVGFR(sel_ix1,psth,bw,toi)
FR = cell(1,length(sel_ix1));

dt = [toi(1):bw:toi(2)];

for it = 1:length(sel_ix1)
    
    T = dt(2)-dt(1);
    tInt = [toi(1):bw:toi(2)];
    x = [];
    for ft = 1:length(tInt)-1
        nSpk = sum(psth{it}(:,find(dt >=tInt(ft) & dt < tInt(ft+1))),2);
        x(:,ft) = nSpk./T;% poststim spike count for each trial, note this is not firing rate
    end;
    FR{it} = mean(x,2);
           
end;
