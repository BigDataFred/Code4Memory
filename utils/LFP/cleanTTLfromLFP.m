function [dataSamplesClean] = cleanARTEFACTfromLFP(dataSamples,ARTix,dt,step,model, Fs)
%
% [Inputs]: dataSamples = your signal (eg LFP, MEG, EEG, etc), must have dim 1xn
%         ARTix       = sample indexes at which artefact occurs
%         dt          = width of the epoch in ms over which the artefact occurs, 
%                       note that units must be in samples (eg 1ms = 32 samples if Fs = 32kHZ) 
%         step        = number of events that are counted for computing the
%                       running template (eg n of events included in
%                       average)
%        model        = should be set to 1, ( use 2 for a high-pass filtered version )
%        Fs           = the sampling frequency of the data (units must be in Hz)
%
% [Output]: dataSamplesClean = cleaned version of the LFP signal.

% F.Roux, UoB
% Birmingham, 2017

if model ==1
    x = [];
    for it =1:length(ARTix)
        x(it,:) = dataSamples(ARTix(it)-dt(1):ARTix(it)+dt(2));
    end;
else
    [b,a] = butter(2,300/(Fs/2),'high');% apply low-pass for LFP
    dum = filtfilt(b,a,dataSamples);
    x2 = [];
    for it = 1:length(ARTix)
        [x2(it,:)] = dum(ARTix(it)-dt(1):ARTix(it)+dt(2));
    end;
end;

ix = 1:step;
c = 0;
dataSamplesClean = dataSamples;
for kt = 1:length(ARTix)
    c = c+1;
    
    Y = dataSamples(ARTix(kt)-dt(1):ARTix(kt)+dt(2));
    Y = Y';
    Y = Y-mean(Y);
    
    if model ==1
        X =mean(x(ix,:),1);
        X = X';
        
        b = regress(Y,[ones(size(X)) X X.^2 X.^3]);
        Yp = [ones(size(X)) X X.^2 X.^3]*b;
    else
        X =mean(x2(ix,:),1);
        X = X';
        
        b = regress(Y,[ones(size(X)) X X.^2 X.^3]);
        Yp = [ones(size(X)) X X.^2 X.^3]*b;
    end;
    
    if c ==step
        if ix(end)+step < length(ARTix)
            ix = ix+step;
        else
            ix = ix(end)+1:length(ARTix);
        end;
        c = 0;
    end;
    
    dataSamplesClean(ARTix(kt)-dt(1):ARTix(kt)+dt(2)) = dataSamples(ARTix(kt)-dt(1):ARTix(kt)+dt(2)) - Yp';
    
end;