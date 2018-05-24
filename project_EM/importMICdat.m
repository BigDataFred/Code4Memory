function [micDat] = importMICdat(p2micDat,micFiles)

[micDat] = cell(1,length(micFiles));
for it = 1:length( micFiles )
    fprintf([num2str(it),'/',num2str(length(micFiles))]);
    [dat] = load([p2micDat,micFiles(it).name]);
    micDat{it} = dat;
    fprintf('\n');
end;