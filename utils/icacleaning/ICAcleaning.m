function ICAcleaning(subj_idx)
%%
restoredefaultpath;
addpath('/bcbl/home/home_a-f/froux/fieldtrip-20150801/');
ft_defaults;
%%
if nargin == 0
    subj_idx = 1;%56:66;
end;
%%
path2files = '/bcbl/home/home_a-f/froux/MEG/Pilot12_VSLmeg_adapted2/140609/';
files1 = dir([path2files,'Pilot12_preproc_preICA_ALL_ds.mat']);
files2 = dir([path2files,'Pilot12_b2-b4_ICs.mat']);
%%
for ot = subj_idx%loop over files
    tic;
    load([path2files,files1(ot).name]);%load the raw meg data
    load([path2files,files2(ot).name]);%load the ICA data
    %% computes the spectrum of all ICs
    cfg = [];
    cfg.method = 'mtmfft';
    cfg.taper = 'dpss';
    cfg.tapsmofrq = 1.5;
    cfg.pad = 'maxperlen';
    cfg.foi = 1:99.99;
    cfg.keeptrials = 'yes';
    cfg.output = 'fourier';
    
    [comp_cf1] = ft_freqanalysis(cfg,comp1);
    [comp_cf2] = ft_freqanalysis(cfg,comp2);
    
    cfg = [];
    cfg.keeptrials = 'yes';
    [comp_pow1] = ft_freqdescriptives(cfg,comp_cf1);
    [comp_pow2] = ft_freqdescriptives(cfg,comp_cf2);
    %% do a first pass to detect EKG-components in the ICA data
    cfg =[];
    cfg.ecg_freq_min = 30;
    cfg.ecg_freq_max = 210;
    cfg.treshold = 2;
    
    [ecg_idx1] = ecg_comp_detect2(cfg,comp1);
    [ecg_idx2] = ecg_comp_detect2(cfg,comp2);
    %% do second pass to detect EKG-components in the ICA data
    cfg = [];
    cfg.ecg_idx = ecg_idx1;
    cfg.meg_data = meg_data1;
    
    [ecg_idx1,qrsA] = ft_component_class_ekg2(cfg,comp1);
    
    cfg = [];
    cfg.ecg_idx = ecg_idx2;
    cfg.meg_data = meg_data2;
    
    [ecg_idx2,qrsA] = ft_component_class_ekg2(cfg,comp2);
    %% do third pass to detect EKG-components in the ICA data
    cfg = [];
    cfg.ecg_idx = ecg_idx1;
    [ecg_idxa] = ecg_plv_detect(cfg,comp_cf1);    
    
    cfg = [];
    cfg.ecg_idx = ecg_idx2;
    [ecg_idxb] = ecg_plv_detect(cfg,comp_cf2);  
    %% do fourth pass to detect EKG-components in the ICA data
    if ~isempty(ecg_idx2)
        cfg = [];
        cfg.ecg_idx = ecg_idx2(find(~ismember(ecg_idx2,ecg_idx)));
        cfg.meg_data = meg_data;
    
        [ecg_idx2,qrsA2] = ft_component_class_ekg2(cfg,comp);
        
        cfg =[];
        cfg.ecg_freq_min = 30;
        cfg.ecg_freq_max = 210;
        cfg.treshold = 2;
        cfg.ecg_indx = ecg_idx2;
        
        [ecg_idx2] = ecg_comp_detect2(cfg,comp);
        
        if size(ecg_idx2,2) > size(ecg_idx2,1)
            ecg_idx2 = ecg_idx2';
        end;
    end;
    %%    
    if ~isempty(ecg_idx2)
        ecg_idx = [ecg_idx;ecg_idx2];
        qrsA(end+1:end+length(qrsA2)) = qrsA2;
        [ecg_idx,s_idx] = sort(ecg_idx);
        qrsA = qrsA(s_idx);
    end;   
    %% visualize ICs carrying ECG activity
    %visualize_ecg_ICs(meg_data,comp,ecg_idx,qrsA,comp_pow);
    %saveas(gcf,['/data/home1/froux/DATASETS/DevelopmentData/figures/ICA_cleaning/',files1(ot).name(1:5),'_ECG_comp.fig'],'fig');
    %close all;
    %%          
    [eog_idx] = eog_comp_detect2(comp,eeg_data);
    %%
    cfg = [];
    cfg.eog_idx = eog_idx;
    [eog_idx2,bData] = eog_comp_detect3(cfg,comp);
  
    if size(eog_idx2,2)>size(eog_idx2,1)
        eog_idx2 = eog_idx2';
    end;
    eog_idx = sort([eog_idx;eog_idx2]);
    %visualize_eog_ICs(eog_idx,comp,comp_pow,bData);
   %saveas(gcf,['/data/home1/froux/DATASETS/DevelopmentData/figures/ICA_cleaning/',files1(ot).name(1:5),'_EOG_comp.fig'],'fig');
   %close all;
    %%
    [SSpidx,protoBlink,bcComp] = SSp_comp_detect(eeg_data,comp);
    %visualize_SSp_ICs(SSpidx,protoBlink,bcComp,comp);
    %saveas(gcf,['/data/home1/froux/DATASETS/DevelopmentData/figures/ICA_cleaning/',files1(ot).name(1:5),'_SSp_comp.fig'],'fig');
    %close all;
    %%
    [emg_idx] = emg_comp_detect(comp_pow);    
    %%
    cfg = [];
    cfg.frequency = [30 100];
    
    comp_freq2 = ft_selectdata(cfg,comp_cf);

    [emg_idx2] = emg_coh_detect(emg_idx,comp_freq2);
    %%
    [emg_idx] = sort([emg_idx';emg_idx2]);    
    %%
    cfg = [];
    cfg.keeptrials = 'no';
    
    [comp_pow] = ft_freqdescriptives(cfg,comp_cf);
    
    %visualize_emg_ICs(comp_pow,comp,emg_idx);
    %saveas(gcf,['/data/home1/froux/DATASETS/DevelopmentData/figures/ICA_cleaning/',files1(ot).name(1:5),'_EMG_comp.fig'],'fig');
    %close all;
    %%
    [line_idx] = line_freq_comp_detect(comp_pow);
    %visualize_line_freq_comp(comp,comp_pow,line_idx);
    %saveas(gcf,['/data/home1/froux/DATASETS/DevelopmentData/figures/ICA_cleaning/',files1(ot).name(1:5),'_linefreq_comp.fig'],'fig');  
    %close all;
    %%
    if size(ecg_idx,2) > size(ecg_idx,1)
        ecg_idx = ecg_idx';
    end;
    if size(eog_idx,2) > size(eog_idx,1)
        eog_idx = eog_idx';
    end;
    if size(SSpidx,2) > size(SSpidx,1)
        SSpidx = SSpidx';
    end;
    if size(emg_idx,2) > size(emg_idx,1)
        emg_idx = emg_idx';
    end;
    if size(line_idx,2) > size(line_idx,1)
        line_idx = line_idx';
    end;
    save([path2files,files1(ot).name(1:5),'_ICartefacts_IDX.mat'],'ecg_idx','eog_idx','emg_idx','SSpidx','line_idx');
    
    comp_idx = unique(sort([ecg_idx;eog_idx;SSpidx;emg_idx;line_idx]));
    cfg = [];
    cfg.component = unique(sort(comp_idx));

    meg_data = ft_rejectcomponent(cfg,comp,meg_data);

    save([path2files,files1(ot).name(1:end-15),'postICA_',date,'.mat'],'meg_data');
    toc;
end;
exit;