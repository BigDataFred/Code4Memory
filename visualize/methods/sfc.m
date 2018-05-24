%%
clear;
clc;

Fs = 1000;
t = 0:1/Fs:6;
f = 20;

lfp= [];
spk = [];
phi = [];
for it = 1:100
    lfp(:,it) = (0+.35*randn(1,length(t)));%sin(2*pi*10.*t)+
    lfp(:,it) = lfp(:,it) - mean(lfp(:,it));
    
    phi(:,it) = angle(hilbert(lfp(:,it)));
    
    [~,ix] = localmax(lfp(:,it)');
    ix = ix(sign(lfp(ix,it))==1);
    
    n = length(ix);
    pct = ceil(n/100*50);
    sel = randperm(n);
    sel = sel(1:pct);
    ix = ix(sel);
    
    
    if it ==1
        ix2 = ix;
    end;
    
    spk(:,it) = zeros(1,length(t));
    spk(ix,it) = 1;
end;
lfp = lfp+spk;

%lfp_lo = lfp;
fNq = Fs/2;
Wn = 300./fNq;
ord = 100;
b = fir1(ord,Wn);

N = size(lfp,1);
lfp_lo = zeros(size(lfp));
for it = 1:size(lfp,2)
    x = lfp(:,it);
    x = [fliplr(x) x fliplr(x)];
    x = filtfilt(b,1,x);
    x = x(N+1:N*2);
    lfp_lo(:,it) = x;
    
end;

dt = [2 8];
N = size(lfp,1);
lfp_lo2 = lfp_lo;
for it = 1:size(lfp,2)
    
    sel = find(spk(:,it));   
    for kt = 1:length(sel)
        spk_ix = sel(kt)-dt(1):sel(kt)+dt(2);
        spk_ix = spk_ix(spk_ix >(dt(1)+2) & spk_ix <=N-(dt(2)+2));               
        
        if ~isempty(spk_ix)
            sel2 = [spk_ix(1)-2:spk_ix(1)-1 spk_ix(end)+1:spk_ix(end)+2];
            
            y = lfp_lo(sel2,it);
            x = sel2;
        
%             y(spk_ix) = [];
%             x(spk_ix) = [];
            
            yi = interp1(x,y,spk_ix,'linear');
            if any(isnan(yi))
                error('NaN detected');
            end;
            lfp_lo2(spk_ix,it) = yi;            
        end;
        
    end;                
    
end;

win = 1e3;
STA = zeros(win*2+1,size(spk,2));
STA2 = zeros(win*2+1,size(spk,2));
c = 0;
for jt = 1:size(spk,2)
    ix = find(spk(:,jt)~=0);
    for kt = 1:length(ix)
        if ((ix(kt))-win >0) && ((ix(kt))+win <=size(spk,1))
            STA(:,jt) = STA(:,jt) + lfp(ix(kt)-win:ix(kt)+win,jt)./length(ix);
            STA2(:,jt) = STA2(:,jt) + lfp_lo2(ix(kt)-win:ix(kt)+win,jt)./length(ix);
        end;
    end;
end;


figure;
subplot(311);
hold on;
plot(t,lfp_lo2(:,1),'r','LineWidth',3);
plot(t,lfp(:,1));
plot(t(ix2),lfp(ix2,1),'ro');
xlim([3 3.4]);
subplot(312);
hold on;
plot(t,phi(:,1),'b.');
plot(t(ix2),phi(ix2,1),'ro');
xlim([3 3.4]);
subplot(313);
hold on;
plot(-win:win,mean(STA,2));
plot(-win:win,mean(STA2,2),'r')
axis tight;xlim([-150 150]);

%%
TW = 8;
k = 2*TW-1;

params                  = [];
params.Fs               = Fs;
params.pad              = 0;
params.fpass            = [1 500];
params.tapers           = [TW k];
params.err              = [2 0.05];
params.trialave         = 1;

[C,~,~,Syy,Snn,f2]=coherencycpb(lfp_lo,spk,params,0);

movingwin = [.5 .025];

TW = 4;
k = 2*TW-1;

params                  = [];
params.Fs               = Fs;
params.pad              = 0;
params.fpass            = [1 500];
params.tapers           = [TW k];
params.err              = [2 0.05];
params.trialave         = 1;
      
[Cgrm,~,~,~,~,tAx,fAx]=cohgramcpb(lfp_lo,spk,movingwin,params,0);

figure;
subplot(2,2,1);
plot(f2,Syy);
axis tight;

subplot(2,2,2);
plot(f2,C);
axis tight;

subplot(2,2,3:4);
imagesc(tAx,fAx,Cgrm');
axis xy;
xlim([3 3.4]);

%%
TW = 8;
k = 2*TW-1;

params                  = [];
params.Fs               = Fs;
params.pad              = 0;
params.fpass            = [1 500];
params.tapers           = [TW k];
params.err              = [2 0.05];
params.trialave         = 1;

[C,~,~,Syy,Snn,f2]=coherencycpb(lfp_lo2,spk,params,0);

movingwin = [.5 .025];

TW = 4;
k = 2*TW-1;

params                  = [];
params.Fs               = Fs;
params.pad              = 0;
params.fpass            = [1 500];
params.tapers           = [TW k];
params.err              = [2 0.05];
params.trialave         = 1;
      
[Cgrm,~,~,~,~,tAx,fAx]=cohgramcpb(lfp_lo2,spk,movingwin,params,0);

figure;
subplot(2,2,1);
plot(f2,Syy);
axis tight;

subplot(2,2,2);
plot(f2,C);
axis tight;

subplot(2,2,3:4);
imagesc(tAx,fAx,Cgrm');
axis xy;
xlim([3 3.4]);
%%
Fs = 1000;
t = 0:1/Fs:6;

lfp = 0+.1*randn(length(t),1);

spk = zeros(length(t),1);
g = 4;
z = [0 1];
for it = 1:length(t)
    
    f = randperm(2);
    
    if z(f(1)) ==1
        if it >3 && sum(spk(it-3:it-1))~=0
        else
            spk(it) = 1;
        end;
    end;
    
end;

ts = t(find(spk)).*1e3;

[isi_count] = hist(diff(ts),1:250);

figure;
plot(1:250,isi_count);

TW = 8;
k = 2*TW-1;

params                  = [];
params.Fs               = Fs;
params.pad              = 0;
params.fpass            = [1 500];
params.tapers           = [TW k];
params.err              = [2 0.05];
params.trialave         = 1;

[C,~,~,Syy,Snn,f2]=coherencycpb(lfp,spk,params,0);

movingwin = [.5 .025];

TW = 4;
k = 2*TW-1;

params                  = [];
params.Fs               = Fs;
params.pad              = 0;
params.fpass            = [1 500];
params.tapers           = [TW k];
params.err              = [2 0.05];
params.trialave         = 1;
      
[Cgrm,~,~,~,~,tAx,fAx]=cohgramcpb(lfp,spk,movingwin,params,0);

figure;
subplot(2,2,1);
plot(f2,Syy);
axis tight;

subplot(2,2,2);
plot(f2,C);
axis tight;

subplot(2,2,3:4);
imagesc(tAx,fAx,Cgrm');
axis xy;