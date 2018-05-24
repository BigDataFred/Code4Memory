function [cleanSignal,noise] = CleanLineNoise(signal,varargin)

%noiseFreq in Hz; 
%Fs (sampling frequency) in s (default Fs = 1000 samples/s;)
%windowSize in s. (default windowSize = 15s)

Fs = 32e3;          %sampling frequency, in Hz. Default 1000
noiseFreq = 50;     %pickup freq to be removed. Default: line noise, 60 Hz
windowSize = 15;    %moving window in sec. Default: 15 sec

for n = 1:2:length(varargin)
    switch varargin{n}
        case 'Fs'
            Fs = varargin{n+1};
        case 'noiseFreq'
            noiseFreq = varargin{n+1};
        case 'windowSize'
            windowSize = varargin{n+1};
        otherwise
            error([mfilename ': Unknown parameter ' varargin{n}]);
    end
end

signalLength = size(signal,2);
signal_timeVector = 1:signalLength;%(1000/Fs):((signalLength-1)*(1000/Fs)+1);

phase = mod(signal_timeVector*noiseFreq/Fs,1.0);
phase = round(phase*32)/32;
[u_phase,~,ind_phase] = unique(phase);
noise = zeros(1,signalLength);

kernel = zeros(size(u_phase));
N = kernel;
windowSize_half = (windowSize*Fs/2);

for p = 1:signalLength
    if p==1
        for o = 1:length(u_phase)
            ind_w = 1:(windowSize_half*2);
            N(o) = sum(ind_phase(ind_w)==o);
            kernel(o) = sum(signal(ind_phase(ind_w)==o));
        end
    elseif ((signalLength-p)>windowSize_half) && (p>windowSize_half) % move kernel
        oldPoint = p-windowSize_half;
        newPoint = p+windowSize_half;
        kernel(ind_phase(oldPoint)) = kernel(ind_phase(oldPoint)) - signal(oldPoint);
        N(ind_phase(oldPoint)) = N(ind_phase(oldPoint)) - 1;
        
        kernel(ind_phase(newPoint)) = kernel(ind_phase(newPoint)) + signal(newPoint);
        N(ind_phase(newPoint)) = N(ind_phase(newPoint)) + 1;
    end
    %gets the average lfp voltage
    noise(p) = kernel(ind_phase(p))./N(ind_phase(p));
    
end
cleanSignal = signal-noise;

% %% sanity check
% nfft = 2^nextpow2(length(noise));
% PSD = fft(noise,nfft)/Fs;
% PSD = abs(PSD).^2;
% 
% PSD = fftshift(PSD);
% PSD = PSD(nfft/2:end);
% 
% nfft = 2^nextpow2(length(signal));
% PSD2 = fft(signal,nfft)/Fs;
% PSD2 = abs(PSD2).^2;
% 
% PSD2 = fftshift(PSD2);
% PSD2 = PSD2(nfft/2:end);
% 
% nfft = 2^nextpow2(length(cleanSignal));
% PSD3 = fft(cleanSignal,nfft)/Fs;
% PSD3 = abs(PSD3).^2;
% 
% PSD3 = fftshift(PSD3);
% PSD3 = PSD3(nfft/2:end);
% 
% f = Fs/2*linspace(0,1,nfft/2+1);
% 
% figure;
% hold on;
% plot(f,log10(PSD));%noise
% plot(f,log10(PSD2),'r');%signal
% plot(f,log10(PSD3),'k');%noise - signal

end
