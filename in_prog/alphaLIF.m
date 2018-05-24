%%
x = zeros(1,200);
x([15 17 19 25 30 45])=1;
T = 20;
tf = find(x==1);
e = zeros(1,length(x));
for kt = 1:length(tf)
for t = 1:length(x)
    
    if t<tf(kt)
    else
        e(t) = e(t)+exp(-(t-tf(kt))/T);
    end;
end;
end;

figure;
subplot(211);
stem(x);
subplot(212);
plot(e);