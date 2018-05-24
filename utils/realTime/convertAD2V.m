function [V] = convertAD2V(dataArray,timeStampArray,numRecordsReturned)
%%
b = length(dataArray)/numRecordsReturned;
dat = reshape(dataArray,[b numRecordsReturned]);
%%
ADBitVolts = 7.6296e-08;

[v,indx] = sort(timeStampArray);
[A,I] = unique(v);
indx = indx(I);

V = double(dat(:,indx)).*ADBitVolts*1e6;