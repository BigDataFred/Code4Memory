function [aIx,dIx] = interpPhaseQuadrants(ft,ix1,ix3)
% estimate ascending and descending phase
ix2 = [];
ix4 = [];
for kt = 1:length(ix1)
    if any(sign(ft(ix1(kt):end))==1)
        ix2(kt) = ix1(kt)+min(find(sign(ft(ix1(kt)+1:end))==1));% ascending zero crossing
    else
        ix2(kt) = NaN;
    end;
    if any(sign(ft(ix3(kt):end))==-1)
        ix4(kt) = ix3(kt)+min(find(sign(ft(ix3(kt)+1:end))==-1));% descending zero crossing
    else
        ix4(kt) = NaN;
    end;
end;
aIx = ix2;
dIx = ix4;