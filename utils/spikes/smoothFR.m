function [s] = smoothFR(y,db,wf)

k = eval([wf,'(db)']);
%k = k./sum(k);
%k = k - mean(k);

s = conv(y,k,'same');

return;
