function [BFlab,BFid,hemiLab,hemiID] = extractMWLabel(chanLab)

ix = regexp(chanLab,'\d{1}');
ix = [ix{:}]-2;
[BFlab] = cell(length(ix),1);
[hemiLab] = cell(length(ix),1);

for it = 1:length(ix)
    BFlab(it) = { chanLab{it}(1:ix(it)) };
    hemiLab(it) = {BFlab{it}(end)};
end;
[BFid] = unique(BFlab);
[hemiID] = unique(hemiLab);