function [Jc] = preData4AversenTestc(data,params)
%%
%TW = 5;
tapers                  = params.tapers;%[TW 2*TW-1];
pad                     = params.pad;%2;
Fs                      = params.Fs;%1e3;
fpass                   = params.fpass;%[0 30];


data=change_row_to_column(data);
N=size(data,1);
nfft=max(2^(nextpow2(N)+pad),N);
[f,findx]=getfgrid(Fs,nfft,fpass); 
tapers=dpsschk(tapers,N,Fs); % check tapers
Jc=mtfftc(data,tapers,nfft,Fs);
Jc=Jc(findx,:,:);

return;