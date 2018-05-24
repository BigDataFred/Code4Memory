function [ixON,ixOFF] = detectIED(sig,filtsig,thr,ddt,plt)

%% create toy data for testing
if nargin ==0
    Fs = 1e3;
    t = 0:1/Fs:1;
    sig = .1*randn(1,length(t));
    ix = randperm(length(sig));
    sig(ix(1:10)) = 1; 
    thr = 3.5;
end;

%% zscore 

zsc = filtsig;%zscore(filtsig);
thrZ = median(zsc)+thr*iqr(zsc);

[thSig1] = abs(zsc) > thrZ;
% 

%% estimate slope

dt = gradient(filtsig);
thrS = median(dt)+thr*iqr(dt);

[thSig2] = abs(dt) > thrS;

% %% estimate kurtosis
% 
% [kts] = zeros( 1 , length(filtsig) );
% parfor it = 1:length(filtsig)
%     
%     if (it-175>0) && (it+175 < length(filtsig))
%         winIx = it-175:it+175;
%     elseif (it-175 >0)
%         winIx = it-175:it;
%     else
%         winIx = it:it+175;
%     end;
%     
%     kts(it) = kurtosis(abs(filtsig(winIx)));
% 
% end;
% 
% thrS = median(kts)+thr*iqr(kts);
% [thSig3] = kts > thrS./2;

%% measure duration

thSig1 = conv(single(thSig1),boxcar(10),'same');
thSig2 = conv(single(thSig2),boxcar(10),'same');
%thSig3 = conv(single(thSig3),boxcar(10),'same');

thSig1(thSig1~=0) = thSig1(thSig1~=0)./thSig1(thSig1~=0);
thSig2(thSig2~=0) = thSig2(thSig2~=0)./thSig2(thSig2~=0);
%thSig3(thSig3~=0) = thSig3(thSig3~=0)./thSig3(thSig3~=0);

sigIx = thSig1 +thSig2;%thSig1 + thSig2 + thSig3;

sigIx = single(sigIx ==1);
sigIx = conv(sigIx,gausswin(10),'same')./sum(gausswin(10));

%%
ixON = find(sign(gradient(sigIx))==1);
ixOFF = find(sign(gradient(sigIx))==-1);
ixON = fliplr(ixON);
ixON(find(diff(ixON) ==-1)) = [];
ixOFF(diff(ixOFF) ==1) = [];
ixON = fliplr(ixON);

ixON2 = ixON;
ixOFF2 = ixOFF;
ixON2(ixON>ddt) = ixON2(ixON>ddt) - ddt;
ixOFF2(ixOFF<length(sigIx)-ddt*2.5) = ixOFF2(ixOFF<length(sigIx)-ddt*2.5)+ddt*2.5;

%%
dum = sig;
for it = 1:length(ixON)    
    sig(ixON2(it):ixOFF2(it)) = NaN;    
end;

%%
% n = length(zsc);
% if strcmp(plt,'y')
%     subplot(5,1,1);
%     hold on;
%     plot(dum);
%     for kt = 1:length(ixON)
%         ix = ixON2(kt):ixOFF2(kt);
%         plot(ix,dum(ix),'r');
%     end;
%     axis tight;
%     subplot(5,1,2);
%     [ax,~,~] = plotyy(1:n,thSig1,1:n,abs(zsc));
%     axis(ax,'tight');
%     set(ax(1),'YLim',[-.05 1.1]);
%     subplot(5,1,3);
%     [ax,~,~] = plotyy(1:n,thSig2,1:n,abs(dt));
%     axis(ax,'tight');
%     set(ax(1),'YLim',[-.05 1.1]);
%     subplot(5,1,4);
%     [ax,~,~] = plotyy(1:n,thSig3,1:n,kts);
%     axis(ax,'tight');
%     set(ax(1),'YLim',[-.05 1.1]);
%     subplot(5,1,5);
%     hold on;
%     plot(sigIx);    
%     plot(ixON,sigIx(ixON),'r^');
%     plot(ixOFF,sigIx(ixOFF),'b^');
%     axis tight;
% end;


% d = zeros(1,length(ixON));
% for it = 1:length(ixON)
%     try
%     d(it) = abs(ixON(it)-ixOFF(it))/Fs;
%     catch
%         return
%     end;
% end;
% 
% %%
% pow = cell(1, length( ixON ));
% for it = 1:length( ixON )
%     
%     e = sig(ixON(it)-1:ixOFF(it)+1);
%     nfft = 1024;
%     y = fft(e,nfft)./Fs;
%     y = fftshift(y);
%     y = y.*conj(y);
%     y = y(nfft/2:length(y));
%     
%     pow{it} = y;
% end;
% f = Fs/2*linspace(0,1,nfft/2+1);

% %%
% figure;
% for it = 1:length( ixON )
%     
%     e = sig(ixON(it)-1:ixOFF(it)+1);
%     subplot(ceil(length(ixON)/5),5,it);
%     plot(e);
%     axis tight;
%     
% end;
% 
% figure;
% for it = 1:length( pow )
%     subplot(ceil(length(pow)/5),5,it);
%     plot(f,pow{it});
%     axis tight;
%     
% end;
% 
% %%
% figure;
% subplot(311);
% hold on;
% plot(t,sig);
% plot(t(ixON),sig(ixON),'r*');
% plot(t(ixOFF),sig(ixOFF),'g*');
% subplot(312);
% hold on;
% plot(t,thSig1,'r');
% subplot(313);
% hold on;
% plot(t,thSig2,'g');