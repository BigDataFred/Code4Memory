function preData4AversenTestc(data1,data2)
%%
TW = 5;
tapers = [TW 2*TW-1];
pad = 2;
Fs = 1e3;
fpass = [0 30];


data1=change_row_to_column(data1);
N=size(data1,1);
nfft=max(2^(nextpow2(N)+pad),N);
[f,findx]=getfgrid(Fs,nfft,fpass); 
tapers=dpsschk(tapers,N,Fs); % check tapers
J1=mtfftc(data1,tapers,nfft,Fs);
J1=J1(findx,:,:);

data2=change_row_to_column(data2); % changes data2 stored as a row vector to a column vector
[N,C]=size(data2); % size of data2
K=size(tapers,2); % size of tapers
tapers=tapers(:,:,ones(1,C)); % add channel indices to tapers
H=fft(tapers,nfft,1); % fourier transform of the tapers
Nsp=sum(data2,1); % number of spikes in each channel
Msp=Nsp'./N; % mean rate for each channel
meansp=Msp(:,ones(1,K),ones(1,size(H,1)));  % add taper and frequency indices to meansp
meansp=permute(meansp,[3,2,1]); % permute to get meansp with the same dimensions as H
data2=data2(:,:,ones(1,K));% add taper indices to the data2
data2=permute(data2,[1 3 2]); % permute data2 to be of the same dimensions as H
data2_proj=data2.*tapers; % multiply data2 by the tapers
J2=fft(data2_proj,nfft,1); % fft of projected data2
J2=J2-H.*meansp; % subtract the dc
J2=J2(findx,:,:);

return;