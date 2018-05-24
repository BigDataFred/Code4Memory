function [cleanSpectrum,freqBinz] = oneOverfRegression(freqBinz,Spectrum)

[n1,m1] = size(freqBinz);
[n2,m2] = size(Spectrum);

Spectrum = Spectrum(freqBinz >3);
freqBinz = freqBinz(freqBinz >3);

if n1<m1;freqBinz = freqBinz';end;
if n2 < m2;Spectrum = Spectrum';end;

X = [ones(length(freqBinz),1) log10(freqBinz)];
Y = log10(Spectrum);

b = (X)\Y;
yp = X*b;

[cleanSpectrum] = 10.^(Y-yp);



return;