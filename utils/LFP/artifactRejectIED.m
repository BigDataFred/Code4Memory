function [chanIx,trlIx,ixON,ixOFF] = artifactRejectIED(dat1,dat2)

%%
[chanIx,trlIx,ixON,ixOFF] = removeIEDtrl(dat1,dat2,.25,150,'n');

%%
function [chanIx,trlIx,ixON,ixOFF] = removeIEDtrl(dat1,dat2,pctTrsh,dt,plt)

ntrl = length(dat1.trial);

[dat1,ixON,ixOFF] = runIEDdetection(dat1,dat2,dt,plt);

[pct,trlChck] = calculateIEDPct(ntrl,ixON,pctTrsh,plt);

delIx = find( pct >= pctTrsh );
[chanIx] = setdiff(1:length(dat1.label),delIx);

if isempty(chanIx);chanIx = delIx;end;

[trlIx] = trlChck;


