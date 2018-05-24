function [Jpb] = preData4AversenTestpb(data,params)
%%
%TW = 5;
tapers              = params.tapers;%[TW 2*TW-1];
pad                 = params.pad;%2;
Fs                  = params.Fs;%1e3;
fpass               = params.fpass;%[0 30];

data=change_row_to_column(data);
N=size(data,1);
nfft=max(2^(nextpow2(N)+pad),N);
[f,findx]=getfgrid(Fs,nfft,fpass); 
tapers=dpsschk(tapers,N,Fs); % check tapers

data=change_row_to_column(data); % changes data stored as a row vector to a column vector
[N,C]=size(data); % size of data
K=size(tapers,2); % size of tapers
tapers=tapers(:,:,ones(1,C)); % add channel indices to tapers
H=fft(tapers,nfft,1); % fourier transform of the tapers
Nsp=sum(data,1); % number of spikes in each channel
Msp=Nsp'./N; % mean rate for each channel
meansp=Msp(:,ones(1,K),ones(1,size(H,1)));  % add taper and frequency indices to meansp
meansp=permute(meansp,[3,2,1]); % permute to get meansp with the same dimensions as H
data=data(:,:,ones(1,K));% add taper indices to the data
data=permute(data,[1 3 2]); % permute data to be of the same dimensions as H
data_proj=data.*tapers; % multiply data by the tapers
Jpb=fft(data_proj,nfft,1); % fft of projected data
Jpb=Jpb-H.*meansp; % subtract the dc
Jpb=Jpb(findx,:,:);

return;