%% set the path def
restoredefaultpath;
addpath(genpath('/media/rouxf/rds-share/Fred/code/mcode/custom/'));% custom code
addpath(genpath('/media/rouxf/rds-share/Common/releaseDec2015/'));%needed to read Nlx data
addpath(genpath('/media/rouxf/rds-share/Common/osort-v3-rel/'));%needed for spike detection & sorting
addpath(genpath('/media/rouxf/rds-share/Common/chronux_2_11/'));% needed for spectral analysis
addpath(genpath('/media/rouxf/rds-share/Common/wave_clus/'));%needed for spike detection & sorting
fn = dir('/media/rouxf/rds-share/Common/fieldtrip-*');
addpath(['/media/rouxf/rds-share/Common/',fn.name]);% spectral analysis, signal processing, spike detection
ft_defaults;

%%
p2NLX = '/media/rouxf/rds-share/Archive/MICRO/P07/fVSpEM';
NLXsesh = '2017-05-06_19-13-42';
makeSpikes  = 1;
makeLFP     = 1;
makeMUA     = 0;

%%
CSCchans = dir([p2NLX,filesep,NLXsesh,filesep,'*.ncs'])
CSCidx = 1:8;

%% open pool of workers
if isempty(gcp('nocreate'))
    %parpool(length(params.CSC_idx),'SpmdEnabled',false);
    parpool(36,'SpmdEnabled',false);
end;

%%
lfp_dat = cell(1,length(CSCidx));
spk_dat = cell(1,length(CSCidx));
parfor it = 1:length(CSCidx)
    it
    [data,readme,chan] = analyze_spontaneous_data(p2NLX,NLXsesh,CSCidx(it),makeSpikes,makeLFP,makeMUA);
    
    [lfp_dat{it}] = data{3};    
    [spk_dat{it}] = data{1};
    
end;
    
%%
lfp_dat = ft_appenddata([],lfp_dat{:});

%%
wt = spk_dat{3}.waveformtime;

for it = 1:length(spk_dat)
    
    for jt = 1:length(spk_dat{it}.waveform)
        
        ts = spk_dat{it}.timestamp{jt};
        if ~isempty(ts)
            ts = double(ts);% convert to double
            
            ts = ts./1e6;%    convert to sec
            ts = ts.*1e3;%    convert to ms
            ts = ts - min(ts);
            
            iti = diff(ts);
            
            % bin spike times into 1 ms bins (Fs = 1kHz)
            dt = min(ts):1:max(ts);
            x = hist(ts,dt);
            
            if max(x) ~=1
                error('spike count out of range');
            end;
            
            [xc,lag] = xcorr(x-mean(x),0.5e3,'coeff');
            ix =(length(lag)-1)/2+1;
            
            xc(ix-1:ix+1) = NaN;
            
            p = length(find(iti <3))/length(iti);
            
            wvf = squeeze(spk_dat{it}.waveform{jt});
            
            figure;
            subplot(2,2,1);
            hold on;
            plot(wt,wvf,'Color',[.9 0 0]);
            plot(wt,mean(wvf,2),'Color',[0 0 0]);
            axis tight;
            title(CSCchans(it).name);
            subplot(2,2,2);
            dt = 0:1:0.5e3;
            [n_iti,~] = histc(iti,dt);
            bar(dt,n_iti);
            axis tight;
            subplot(2,2,3);
            bar(lag,xc,'k');
            axis tight;
            set(gca,'XTick',[-5:5].*1e3);
        end;
    end;
    
end;

%%
cfg = [];
cfg.viewmode = 'vertical';

ft_databrowser(cfg,lfp_dat);