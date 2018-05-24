function [selIx] = filterChanels4Label(BFLab,filterLabel)

chck = regexp(BFLab,filterLabel);
cnt = 0;
selIx = [];
for it = 1:length( chck )
    if ~isempty(chck{it})
        cnt = cnt+1;
        selIx(cnt) = it;
    end;
end;
    
return;