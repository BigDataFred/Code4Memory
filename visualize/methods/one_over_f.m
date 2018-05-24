
%%
[b] = regress(10*log10(Syy1),[ones(size(Syy1)) log10(f1')]);
yp = b(1)+b(2)*log10(f1');

figure;
subplot(121);
hold on;
plot(log10(f1),10*log10(Syy1));
plot(log10(f1),yp,'r');
axis tight;

subplot(122);
cSyy1 = 10*log10(Syy1) - yp;
cSyy1 = 10.^cSyy1;
plot(f1,cSyy1);
axis tight;