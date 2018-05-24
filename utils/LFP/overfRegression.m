
function [Snn] = overfRegression(Snn,fx)

X = [ones(size(Snn)) log10(fx)];
Y = Snn;
Y = log10(Y);
b = X\Y;
yp = X*b;

Snn = 10.^(Y-yp);