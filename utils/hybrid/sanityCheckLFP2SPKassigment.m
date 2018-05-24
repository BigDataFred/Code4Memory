function sanityCheckLFP2SPKassigment(MWdat,spkDat,BFinf,data_lfp);

x = {};for it = 1:length(spkDat);x(it) = {spkDat{it}.hdr.label};end;

chck = [BFinf(1,:)' MWdat.label([BFinf{2,:}]) x'];

c1 = sum(strcmp(chck(:,2),chck(:,3))) ;
c2 = sum(strcmp(chck(:,1),chck(:,2))) ;
c3 = sum(strcmp(chck(:,1),chck(:,3))) ;

if any([c1 c2 c3] ~= length(data_lfp))
    error('LFP label assignment is inconsistent');
end;

id = unique(BFinf(1,:));
ix = cell(1,length(id));
n = zeros(1,length(id));
for it = 1:length(id)
    
    ix{it} = find(strcmp(BFinf(1,:),id(it)));
    n(it) = length(ix{it});
end;
ix = ix(n>1);

chck = zeros(1,length(ix));
for it = 1:length(ix)
    
    chck(it) = isequal(data_lfp{ix{it}(1)},data_lfp{ix{it}(2)});
    
end;

if sum(chck) ~= length(ix)
     error('LFP label assignment is inconsistent');
end;

return;
