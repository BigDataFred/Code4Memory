function [dataSamplesClean,res] = cleanARTIFACTfromLFP(dataSamples,ARTix,dt,step,model, Fs)
%
% [Inputs]: dataSamples = your signal (eg LFP, MEG, EEG, etc), must have dim 1xn
%         ARTix       = sample indexes at which artefact occurs
%         dt          = width of the epoch in ms over which the artefact occurs, 
%                       note that units must be in samples (eg 1ms = 32 samples if Fs = 32kHZ) 
%         step        = number of events that are counted for computing the
%                       running template (eg n of events included in
%                       average)
%        model        = should be set to 1, ( use 2 for a high-pass filtered version )
%        Fs           = the sampling frequency of the data (units must be
%                       in Hz), default = [], must be specified if model =2
%
% [Output]: dataSamplesClean = cleaned version of the LFP signal.

% F.Roux, UoB
% Birmingham, 2017

x = [];
selIx = [];
tmp = dataSamples;
c=0;

if model ==2    
    [b,a] = butter(2,[167/(Fs/2) 10e3/(Fs/2)],'bandpass');% apply low-pass for LFP
    tmp = filtfilt(b,a,dataSamples);
end;

for it =1:length(ARTix)
    if (ARTix(it)-dt(1)>1 ) && ( ARTix(it)+dt(2)<length(tmp) )
        c= c+1;
        [x(c,:)] = tmp(ARTix(it)-dt(1):ARTix(it)+dt(2));
        selIx(c) = it;
    end;
end;
clear tmp;
ARTix = ARTix(selIx);


c = 0;
dataSamplesClean = dataSamples;
m = length(ARTix);
res = zeros(length(ARTix),1);
for kt =1:length(ARTix)
    c = c+1;
    
    Y = dataSamples(ARTix(kt)-dt(1):ARTix(kt)+dt(2));
    Y = Y';
    Y = Y-mean(Y);    
    
    if kt <=step
        d1 = (kt-1);
        if 2*step < size(x,1)
            d2 = step*2-d1;
        else
            d2 =size(x,1)-kt;
        end;
        X = x(kt-d1:kt+d2,:);
    elseif kt >= m-step
        d1 = m-kt;
        if 2*step < size(x,1)
            d2 = step*2-d1;
        else
            d2 = size(x,1)-kt;
        end;
        X = x(kt-d2:kt+d1,:);
    else
        X = x( kt-step:kt+step,:);
    end;    
    X = X';
    X = X-(ones(size(X,1),1)*mean(X,1));
    
    pred = mean(X,2);        
    pred = conv(pred,gausswin(11),'same')./sum(gausswin(11));  
    [b] = regress(hilbert(Y),[ones(size(pred,1),1) hilbert(pred)]);% 
    %[b] = regress((Y),[ones(size(pred,1),1) (pred)]);
    Yp = real([ones(size(pred,1),1) pred]*b);%X.^2 X.^3         
      
    [d] = sqrt(sum((X-(Y*ones(1,size(X,2)))).^2,1));
    [~,ix] = min(d);    
    Yp2 = X(:,ix);
    
%     dum = dataSamples;
%     x2 = dum(ARTix(kt)-dt(1):ARTix(kt)+dt(2)); 
%     x2 = x2 - Yp';    
%     dum(ARTix(kt)-dt(1):ARTix(kt)+dt(2)) = x2;
%     x1 = dataSamples(ARTix(kt)-0.5*32:ARTix(kt)+5*32);
%     x2 = dum(ARTix(kt)-0.5*32:ARTix(kt)+5*32);
    
    d1 = sqrt(sum((Y-Yp )).^2);
    
%     dum = dataSamples;
%     x2 = dum(ARTix(kt)-dt(1):ARTix(kt)+dt(2)); 
%     x2 = x2 - Yp2';    
%     dum(ARTix(kt)-dt(1):ARTix(kt)+dt(2)) = x2;
%     x1 = dataSamples(ARTix(kt)-0.5*32:ARTix(kt)+5*32);
%     x2 = dum(ARTix(kt)-0.5*32:ARTix(kt)+5*32);
    
    d2 = sqrt(sum((Y-Yp2)).^2);
    
    if d2< d1
        Yp = Yp2;
    end;
    
    x2 = dataSamplesClean(ARTix(kt)-dt(1):ARTix(kt)+dt(2)); 
    x2 = x2 - Yp';    
    dataSamplesClean(ARTix(kt)-dt(1):ARTix(kt)+dt(2)) = x2;
    
    if ( ARTix(kt)-0.5*32 > 1) && ( ARTix(kt)+5*32 < length(dataSamplesClean) )
        d = dataSamplesClean(ARTix(kt)-0.5*32:ARTix(kt)+5*32).^2;
    elseif ( ARTix(kt)-0.5*32 > 1) && ( ARTix(kt)+5*32 > length(dataSamplesClean) )
        d = dataSamplesClean(ARTix(kt)-0.5*32:end).^2;
    else
        d = dataSamplesClean(1:ARTix(kt)+5*32).^2;
    end;
    
    d = sqrt(sum(d))/length(d);
    res(kt) = d;

%     figure;
%     subplot(131);
%     hold on;
%     plot(X);
%     plot(Y,'r');
%     legend('predictor','signal');
%     subplot(132);
%     hold on;
%     plot(Y,'r');
%     plot(Yp,'g');
%     legend('signal','fit');
%     subplot(133);
%     hold on;
%     plot(dataSamples(ARTix(kt)-dt(1)*2:ARTix(kt)+dt(2)*2));
%     plot(dataSamplesClean(ARTix(kt)-dt(1)*2:ARTix(kt)+dt(2)*2),'r');
%     legend('signal','signal clean');
%     pause;
%     clf;
     
end;