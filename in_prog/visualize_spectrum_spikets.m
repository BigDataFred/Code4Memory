%%
clear all;
clc;
close all;

% %%
% Fs = 1e3;
% x = zeros(1,Fs);
% sel = randperm(Fs);
% sel = sel(1:2);
% x(sel) = 1;
% 
% nfft = 2^nextpow2(Fs);
% y = fft(x,nfft)./Fs;
% y = conj(y).*y;
% y = fftshift(y);
% y = y(nfft/2:end);
% 
% f = Fs/2*linspace(0,1,nfft/2+1);
% 
% figure;
% plot(f,y);

%%
Fs = 1000;
t = 1e-4:1/Fs:1;
dt = 1e-4:1e-3:1;

x2 = sin(2*pi*5.*t);

trl = 100;
X = zeros(length(t),trl);

for it = 1:trl
    
    f = randperm(3);
    f =f(1);
    
    x = zeros(1,length(t));
    
    if f ==1
        
        sel = [];        
        
        sel = find(x2 >= max(x2)/100*90);
        ix = randperm(length(sel));
        ix = ix(1:ceil(length(sel)/1.25));
        sel(ix) = [];
        
        x(sel) = 1;
        
        sel2 = [];
        sel2 = randperm(length(t));
        sel2 = sel2(1:150);
        
        x(sel2) = 1;
                
        [dum] = histc(t(find(x~=0)),dt);
        
        X(:,it) = dum;
        ref = dum;
        
    elseif f ==2
        
        sel2 = [];
        sel2 = randperm(length(t));
        sel2 = sel2(1:150);
        x(sel2) = 1;
                
        [dum] = histc(t(find(x~=0)),dt);
        
        X(:,it) = dum;
    else
        sel2 = [];
        sel2 = randperm(length(t));
        sel2 = sel2(1:150);
        x(sel2) = 1;
                
        [dum] = histc(t(find(x~=0)),dt);
        
        X(:,it) = dum;
        
    end;
    
end;

%figure;
%imagesc(t,1:size(X,2),X');

%%
T = length(x)/Fs;
W =1/T;

params                  = [];
params.tapers           = [W T 1];
params.Fs               = Fs;
params.pad              = 0;
params.fpass            = [1 100];
params.err              = [1 0.05];
params.trialave         = 1;

[S,freq,R,Serr] = mtspectrumpb(X,params,0);

figure;
subplot(121);
imagesc(t,1:size(X,2),X');
%hold on;
%plot(t,ref,'b.');
%plot(t,x2,'r.');
%lim([0 1]);

subplot(122);
plot(freq,S./R,'b');
axis tight;