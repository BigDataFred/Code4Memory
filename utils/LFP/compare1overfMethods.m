%%
[dt] = MWdat.time{1}(2)-MWdat.time{1}(1);
Fs = 1/dt;

movingwin1 = [1 .01];

TW = 3;
params1                  = [];
params1.tapers           = [TW , 2*TW-1];
params1.pad              = -1;
params1.Fs               = Fs;
params1.err              = 0;
params1.trialave         = 1;
params1.fpass            = [0 30];

S  = cell(1,length(data_lfp));
S2  = cell(1,length(data_lfp));
S3  = cell(1,length(data_lfp));

for it = 1:length(data_lfp)
    [S{it},f] = mtspectrumc(data_lfp{it},params1);
    [Fx,Fy] = gradient(data_lfp{it});
    [S2{it}] = mtspectrumc(Fy,params1);
    [S3{it}] = mtspectrumc(data_lfp2{it},params1);
end;
f = f';

[Fx,Fy] = gradient(LFPsig{1});
[S5] = mtspectrumc(Fy,params1);

x = zeros(size(S{1}));
for it = 1:length(S)
    x = x + S{it};
end;
x = x./it;

xII = zeros(size(S2{1}));
xIII = zeros(size(S2{1}));
for it = 1:length(S2)
    xII = xII + S2{it};
    xIII = xIII + S3{it};
end;
xII = xII./it;
xIII = xIII./it;

b = regress(log10(x),[ones(size(x)) log10(f)]);
yp1 = b(1)+b(2).*log10(f);
a1 = b(2);

x2 = x(find(f>=1));
f2 = f(find(f>=1));
b = regress(log10(x2),[ones(size(x2)) log10(f2)]);
yp2 = b(1)+b(2).*log10(f2);
a2 = b(2);

x3 = x(find(f>=2.5));
f3 = f(find(f>=2.5));
b = regress(log10(x3),[ones(size(x3)) log10(f3)]);
yp3 = b(1)+b(2).*log10(f3);
a3 = b(2);

figure;
subplot(231);
hold on;
plot(log10(f),log10(x));
plot(log10(f),yp1,'c');
plot([0 0],[min(yp1) max(yp2)],'r--');
axis tight;
title(['frequency range: ',num2str(round(min(f)*10)/10),'-',num2str(round(max(f)*10)/10)]);
subplot(232);
hold on;
plot(log10(f),log10(x));
plot(log10(f2),log10(x2),'r');
plot(log10(f2),yp2,'g');
plot([0 0],[min(yp1) max(yp2)],'r--');
axis tight;
title(['frequency range: ',num2str(round(min(f2)*10)/10),'-',num2str(round(max(f2)*10)/10)]);
subplot(233);
hold on;
plot(log10(f),log10(x));
plot(log10(f3),log10(x3),'k');
plot(log10(f3),yp3,'m');
plot([0 0],[min(yp1) max(yp2)],'r--');
axis tight;
title(['frequency range: ',num2str(round(min(f3)*10)/10),'-',num2str(round(max(f3)*10)/10)]);

subplot(234);
hold on;
plot(f,x-10.^(yp1));
axis tight;
xlim([0 30]);
title('1/f corrected');
subplot(235);
hold on;
plot(f2,x2-10.^(yp2));
axis tight;
xlim([0 30]);
title('1/f corrected');
subplot(236);
hold on;
plot(f3,x3-10.^(yp3));
axis tight;
xlim([0 max(f)]);
title('1/f corrected');

a = get(gcf,'Children');
for it = 1:length(a);
    xlabel(a(it),'Frequency [Hz]');
    ylabel(a(it),'Power');
end;


figure;
hold on;
plot(f,xII);
plot(f,xIII,'r');
plot(f,S5,'k');
a = get(gcf,'Children');
for it = 1:length(a);
    xlabel(a(it),'Frequency [Hz]');
    ylabel(a(it),'Power');
end;
title('1/f corrected');

%%
TW = 3;
params1                  = [];
params1.tapers           = [TW , 2*TW-1];
params1.pad              = -1;
params1.Fs               = Fs;
params1.err              = 0;
params1.trialave         = 1;
params1.fpass            = [2 30];

SgmLFP  = cell(1,length(data_lfp));
SgmLFP2 = cell(1,length(data_lfp));
SgmLFP3 = cell(1,length(data_lfp));
for it = 1:length(data_lfp)
    fprintf([num2str(it),'/',num2str(length(data_lfp))]);
    [Fx,Fy] = gradient(data_lfp{it});
    [SgmLFP1{it},tx,fx] = mtspecgramc(Fy,movingwin1,params1);
    [SgmLFP2{it},~,~] = mtspecgramc(data_lfp2{it},movingwin1,params1);
    [SgmLFP3{it},~,~] = mtspecgramc(data_lfp{it},movingwin1,params1);
    fprintf('\n');
end;

SgmLFP4 = SgmLFP3;
selIx = 1:length(fx);%find(fx >=1);
%delIx = setdiff(1:length(fx),selIx);
for it = 1:length(SgmLFP3)
    
    lf = fx(selIx);
    lf = log10(lf)';
    X = [ones(size(lf)) lf];
    for jt = 1:size(SgmLFP3{it},1)
        lp = log10(SgmLFP3{it}(jt,selIx))';
        b = X\lp;
        yp = X*b;
        SgmLFP4{it}(jt,selIx) = 10.^(lp-yp)';
        %SgmLFP4{it}(jt,delIx) = 0;
    end;
        
end;

% fx(delIx) = [];
% for it = 1:length(SgmLFP)   
%     SgmLFP{it}(:,delIx) = [];
%     SgmLFP2{it}(:,delIx) = [];
%     SgmLFP3{it}(:,delIx) = [];
% end;

figure;
cnt1 = 0;
cnt2 = 8;
cnt3 = 16;
cnt4 = 24;
for it = 1:length(SgmLFP)
    
    cnt1 = cnt1+1;
    subplot(4,8,cnt1);
    imagesc(tx-2,fx,(SgmLFP3{it})');
    axis xy;
    
    cnt2 = cnt2+1;
    subplot(4,8,cnt2);
    imagesc(tx-2,fx,(SgmLFP1{it})');
    axis xy;
    
    cnt3 = cnt3+1;
    subplot(4,8,cnt3);
    imagesc(tx-2,fx,(SgmLFP2{it})');
    axis xy;
        
    cnt4 = cnt4+1;
    subplot(4,8,cnt4);
    imagesc(tx-2,fx,(SgmLFP4{it})');
    axis xy;
    
end;