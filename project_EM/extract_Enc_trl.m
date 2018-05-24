function [trl_all] = extract_Enc_trl(LogDat)

%find those trials that correspond to Encoding
dum = [];
dum(1,:) = LogDat.idx(1,:);
for it = 2:size(LogDat.idx,1)
    dum(it,:) = [LogDat.idx(it,:)+max(LogDat.idx(it-1,:))]
end;

%sanity check
if diff(LogDat.idx,[],2) ~= diff(dum,[],2)
    error('trial assignment must match');
end;

% generate vector with indices for each run
dum2 = cell(1,size(dum,1));
for it = 1:size(dum,1)
    dum2{it} = dum(it,1):dum(it,2);% from start to end of run
end;

trl_all = [dum2{:}];% all trials