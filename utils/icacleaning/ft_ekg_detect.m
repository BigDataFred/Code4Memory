function ft_ekg_detect(comp)

%% do a first pass to detect EKG-components in the ICA data
[ekg_idx1] = ecg_comp_detect(comp);
%% test phase stability of detected EKG components
k = 0;% counter variable
qrsA = {};% is needed for visualization further downstream in pipeline
artifact = cell(1,length(ekg_idx1));
comp_idx1 = zeros(length(ekg_idx1),1);
qrsA = cell(length(ekg_idx1),1);
dum = struct;
dum.time = comp.time;
dum.trial = comp.trial;
dum.label = meg_data.label;
dum.sampleinfo = meg_data.sampleinfo;
for jt = 1:length(ekg_idx1)
    %this section of code creates a dummy structure with a
    %concatednated signal of the EKG components detected during the
    %first-pass
    
    cfg.trl = zeros(length(dum.trial),3);
    concat_ekg = zeros(1,length(dum.trial)*length(dum.time{1}));
    idx = 1:length(dum.time{1});
    for kt = 1:length(dum.trial)
        cfg.trl(kt,1:2) = [1+length(dum.trial{kt})*(kt-1) length(dum.trial{kt})+length(dum.trial{kt})*(kt-1)];
        cfg.trl(kt,3) = cfg.trl(kt,1)+find(dum.time{kt} == 0);
        concat_ekg(idx) = dum.trial{kt}(ekg_idx1(jt),:);
        idx = idx +length(dum.time{kt});
    end;
    
    dum.sampleinfo = cfg.trl(:,1:2);
    
    %this section searches for QRS-like waveform and returns sample
    %indexes
    cfg = [];
    cfg.trl        = [dum.sampleinfo zeros(size(dum.sampleinfo,1),1)];
    cfg.continuous = 'no';
    cfg.artfctdef.ecg.channel = dum.label{ekg_idx1(jt)};
    cfg.artfctdef.ecg.feedback =  'no';
    
    % store the sample information of QRS events for each detected EKG
    % component
    [cfg, artifact{jt}] = ft_artifact_ecg(cfg, dum);
    
    if (length(artifact{jt}) < 500)
        artifact{jt} = [];
    else
        %% organize QRS-locked ICA data
        QRS = zeros(size(artifact{jt},1),length(artifact{jt}(1,1):artifact{jt}(1,2)));
        hQRS = zeros(size(artifact{jt},1),length(artifact{jt}(1,1):artifact{jt}(1,2)));
        for kt = 1:size(artifact{jt},1)
            QRS(kt,:) = concat_ekg(artifact{jt}(kt,1):artifact{jt}(kt,2));%Amplitude information
            hQRS(kt,:) = hilbert(concat_ekg(artifact{jt}(kt,1):artifact{jt}(kt,2)));% Analytical signal - needed for phase estimation
        end;
        
        % measure the stability of the phase angle over subsequent
        % QRS-waves
        itc = [];
        itc = abs(sum(hQRS,1))./size(hQRS,1);
        
        %measure the amplitude and maximum of the average QRS-locked
        %wave
        A = [];
        A =squeeze(mean((QRS),1));
        pix1 = find(A == max(A));
        
        %measure KLD between the PDF of phase angles and the uniform
        %ditribution
        n1 = [];
        [n1,dum2] = hist(angle(hQRS(:,pix1)),[-pi:pi/6:pi]);
        n1 = n1./sum(n1);
        n1 = n1+1e-12;
        h = -sum(n1.*log(n1),2);
        kl1 = log(length(n1))-h;
        
        %measure the degree of correlation between the average QRS-wave
        %and the single QRS-events
        r = zeros(size(QRS,1),1);
        for nt = 1:size(QRS,1);
            r(nt) = corr(abs(QRS(nt,:))',mean(abs(QRS),1)');
        end;
        p1 = [];
        p1 = 1-length(find(abs(r)>=.5))/length(r);%quantify percentage of stereotypical QRS-events
        
        %measure KLD between the PDF of random phase angles and the uniform
        %ditribution
        kl2 = zeros(200,1);
        shuf_hQRS = zeros(size(QRS,1),size(QRS,2));
        for nt = 1:200
            for mt = 1:size(QRS,1)
                shuf_hQRS(mt,:) = hilbert(QRS(mt,randperm(size(QRS,2))));
            end;
            n2 = [];
            [n2,dum2] = hist(angle(shuf_hQRS(:,pix1)),[-pi:pi/6:pi]);
            n2 = n2./sum(n2);
            n2 = n2+1e-12;
            h = -sum(n2.*log(n2),2);
            kl2(nt) = log(length(n2))-h;
        end;
        p2 = [];
        p2 = length(find(kl2>=kl1))/length(kl2);%quantify distance of empirical phase distribution vs random phase distribution
        
        if p1 < 0.5 && p2 < 0.05%if QRS-wave shows a high degree of phase consistency across trials then update information
            k = k+1;
            % qrsA contains information for downstream visualization of ICs
            % carrying EKG activity
            qrsA{k}.avg = (A-mean(A))./std(A);
            qrsA{k}.trl = QRS;
            qrsA{k}.itc = itc;
            comp_idx1(k) = ekg_idx1(jt);%index refering to ICs carrying EKG activity
        end;
        qrsA(k+1:end) = [];
        comp_idx1(k+1:end) = [];
    end;
end;
