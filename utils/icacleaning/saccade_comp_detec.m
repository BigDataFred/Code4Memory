function [sacc_idx,strAVG,compAVG] = saccade_comp_detec(eeg_data,comp)
%% low-pass filter the H-eog signals
cfg = [];
cfg.channel = eeg_data.label(2);% HEOG 

dum = ft_selectdata(cfg,eeg_data);

cfg = [];
cfg.lpfilter = 'yes';
cfg.lpfreq = 3;

dum = ft_preprocessing(cfg,dum);
%%
win = [zeros(1,dum.fsample*0.65) ones(1,dum.fsample*0.65) zeros(1,dum.fsample*0.65)];%rectangular window
k = 0;
strAVG = zeros(5e3,1.5*dum.fsample+1);
compAVG = zeros(5e3,length(comp.label),1.5*dum.fsample+1);
for it = 1:length(dum.trial)
    
    x = conv(dum.trial{it},win);% convolution with rectangular window
    x(1:length(win)/2) = [];
    x(end-(length(win)/2-2):end) = [];
    
    % z-score convolved signal
    x = x-mean(x);
    x = x./std(x);
    
    % search for trsh-xing
    crsidx = find(x >2);
    crsidx(find(diff(crsidx)==1)+1) = [];        
    
    % compute saccade-triggered average of EEG and ICs data
    for jt = 1:length(crsidx)
        if (crsidx(jt)-dum.fsample*0.75)>0 && (crsidx(jt)+dum.fsample*0.75) < length(dum.trial{it})
            k = k+1;
            strAVG(k,:) = dum.trial{it}(crsidx(jt)-dum.fsample*0.75:crsidx(jt)+dum.fsample*0.75);
            compAVG(k,:,:) = comp.trial{it}(:,crsidx(jt)-dum.fsample*0.75:crsidx(jt)+dum.fsample*0.75);
        end;
    end;
    
end;
strAVG(k+1:end,:) = [];
compAVG(k+1:end,:,:) = [];
%% z-score the saccade triggered Average
strAVG = strAVG - repmat(mean(strAVG,2),[1 size(strAVG,2)]);
strAVG = strAVG ./ repmat(std(strAVG,0,2),[1 size(strAVG,2)]);
%% fake fieldtrip strucutre
dum = struct;
dum.fsample = comp.fsample;
dum.label = comp.label;
dum.time = cell(1,size(compAVG,1));
dum.trial = cell(1,size(compAVG,1));
for it = 1:size(compAVG,1)
    dum.time{it} = -.75:1/dum.fsample:.75;
    dum.trial{it} = squeeze(compAVG(it,:,:));
end;
%% search for ICs that are correlated with H-EOG 
if ~isempty(dum.trial)
    cfg = [];
    cfg.keeptrials = 'no';
    erf1 = ft_timelockanalysis(cfg,dum);
    
    r = zeros(length(erf1.label),1);
    for jt = 1:length(erf1.label)
            [r(jt)] = max(abs(xcorr(erf1.avg(jt,:),mean(strAVG,1)','coeff')));        
    end;
    %%
    sacc_idx = find(abs(r) >.6);
    if isempty(sacc_idx)
        strAVG = [];
        compAVG = [];
    end;
else
    sacc_idx = [];
    strAVG = [];
    compAVG = [];
end;

% %%
% figure;
% subplot(211);
% imagesc(-.75:1/dum.fsample:.75,1:size(strAVG1,1),strAVG2);
% subplot(212);
% plotyy(-.75:1/dum.fsample:.75,mean(strAVG2,1),erf1.time,erf1.avg(find(abs(r) >.65 & p < 1e-9),:));
% %%
% cfg = [];
% cfg.viewmode = 'component';
% cfg.layout = 'neuromag306planar.lay';
% cfg.channel = [13 18];
% 
% ft_databrowser(cfg,comp);