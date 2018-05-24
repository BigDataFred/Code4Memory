function [BFinf] = recodeBFinf(BFinf,origLabel)
%%
ix = zeros(1,length(BFinf(1,:)));
for it = 1:length(BFinf(1,:))
    ix(it) = find(strcmp(BFinf(1,it),origLabel));  
end;

chck = strcmp(BFinf(1,:)',origLabel(ix));

if sum(chck) ~= length(BFinf(1,:))
    error('recoding of channel data is inconsistent');
end;
for it = 1:length(ix)
    BFinf(2,it) = {ix(it)};
end;

return;
