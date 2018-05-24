function [comp_idx,qrs] = ft_component_class_ekg2(cfg,comp)

%% initialize
comp_idx = zeros(length(cfg.ecg_idx),1);
qrs = cell(length(cfg.ecg_idx),1);% is needed for visualization further downstream in pipeline
% counter variables
c1 = 0;
%% concatenate the ICA data across all single trials
nsamples = zeros(length(comp.trial),1);
for kt = 1:length(comp.trial)
    nsamples(kt) = length(comp.trial{kt});
end;
nsamples = sum(nsamples);

concat_ekg = zeros(length(comp.label),nsamples);
idx = 1:length(comp.time{1});
for kt = 1:length(comp.trial)
    concat_ekg(:,idx) = comp.trial{kt};
    if sign(kt-length(comp.trial))==-1
        idx = idx(end)+1:idx(end)+length(comp.trial{kt+1});
    end;
end;
%% this section of code creates a dummy structure with a
%concatednated signal of the EKG components detected during the
%first-pass
dummy2 = struct;
dummy2.time{1} = 0:1/comp.fsample:(size(concat_ekg,2)-1)/comp.fsample;
dummy2.trial{1} = concat_ekg;
dummy2.label = comp.label;
dummy2.fsample = comp.fsample;
dummy2.sampleinfo = [1 length(dummy2.trial{1})];
%% search for QRS-events
for jt = 1:length(cfg.ecg_idx)   
    
    [ekg_sig] = concat_ekg(cfg.ecg_idx(jt),:);
    
    cfg2 = [];
    cfg2.trl        = [dummy2.sampleinfo zeros(size(dummy2.sampleinfo,1),1)];
    cfg2.continuous = 'yes';
    cfg2.artfctdef.ecg.channel = dummy2.label{cfg.ecg_idx(jt)};
    cfg2.artfctdef.ecg.feedback =  'no';
    cfg2.artfctdef.ecg.inspect = dummy2.label{cfg.ecg_idx(jt)};
    
    % store the sample information of QRS events for each detected EKG
    % component    
    [cfg2, artifact] = ft_artifact_ecg(cfg2, dummy2);
       
    sampleinf = zeros(length(comp.trial),2);
    sampleinf(1,:) = [1 length(comp.trial{1})];
    for mt = 2:length(comp.trial)        
        sampleinf(mt,:) = [sampleinf(mt-1,2)+1 sampleinf(mt-1,2)+length(comp.trial{mt})];        
    end;

    ekg_trl = zeros(size(sampleinf,1),1);
    c3 = 0;
    for mt = 1:size(artifact,1)
        for kt = 1:size(sampleinf,1)
            if (artifact(mt,1) >= sampleinf(kt,1)) && (artifact(mt,2) <= sampleinf(kt,2))
                c3 = c3+1;
                ekg_trl(c3) =  kt;
            end;
        end;
    end;
    ekg_trl(c3+1:end) = [];
    
    % compute the % of events per trial duration e.g if the trials was 3 s
    % long and there were 2 events than there were 66% events per second
    trl = unique(ekg_trl);
    n = zeros(length(trl),1);
    for mt = 1:length(trl)
        n(mt) = (length(comp.trial{mt})/comp.fsample)/length(find(ekg_trl==trl(mt)));
    end;
    pct = median(n);
    %% QRS-events need to be present in at least 50 % of all trials 
    if pct >= .65        
        % organize QRS-locked ICA data
        QRS = zeros(size(artifact,1),length(artifact(1,1):artifact(1,2)));
        hQRS = zeros(size(artifact,1),length(artifact(1,1):artifact(1,2)));
        for kt = 1:size(artifact,1)
            QRS(kt,:) = ekg_sig(artifact(kt,1):artifact(kt,2));%Amplitude information
            hQRS(kt,:) = hilbert(ekg_sig(artifact(kt,1):artifact(kt,2)));% Analytical signal - needed for phase estimation
        end;
        
        % measure the stability of the phase angle over subsequent
        % QRS-waves
        itc = [];
        itc = abs(sum(hQRS,1))./size(hQRS,1);
        
        %measure the amplitude and maximum of the average QRS-locked
        %wave
        A = [];
        A =squeeze(mean((abs(QRS)),1));
         
        c1 = c1+1;
        % qrsA contains information for downstream visualization of ICs
        % carrying EKG activity
        qrs{c1}.zavg = (A-mean(A))./std(A);
        qrs{c1}.avg = A;
        qrs{c1}.trl = QRS;
        qrs{c1}.itc = itc;
        qrs{c1}.sampleinfo = artifact;
        
        comp_idx(c1) = cfg.ecg_idx(jt);%index refering to ICs carrying EKG activity

    end;
end;
qrs(c1+1:end) = [];
comp_idx(c1+1:end) = [];
