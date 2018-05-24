function [tw] = createVariableTimeWindowsAcrossBinRange(binRange,increment,maxLength)
%%
%mf = unique(sign(binRange));
%binRange = sort(abs(binRange));
b = binRange(1):increment:binRange(2);
n = length(b)-1;

cnt = 0;
tw = zeros(n*(n-1)/2+n,2);
for it = 1:length(b)-1
    for jt = it+1:length(b)
        cnt = cnt+1;
        tw(cnt,:) = [b(it) b(jt)];
    end;
end;

tw(diff(tw,[],2)<maxLength,:) = [];
% tw = tw.*mf;

return;