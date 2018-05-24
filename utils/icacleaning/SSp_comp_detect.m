function [SSpidx,protoBlink,bcComp] = SSp_comp_detect(eeg_data,comp_data)
%%
SSpidx = [];
protoBlink = [];
bcComp = [];
%% bp-fitler MEG data
cfg = [];
cfg.bpfilt = 'yes';
cfg.bpfreq = [30 (eeg_data.fsample/2)-1];
cfg.bpfilttype = 'but';
cfg.channel = {'all'};

[comp_data] = ft_preprocessing(cfg,comp_data);
%% bp-filter EOG data
cfg = [];
cfg.bpfilt = 'yes';
cfg.bpfreq = [30 (eeg_data.fsample/2)-1];
cfg.bpfilttype = 'but';
cfg.channel = eeg_data.label(1:2);

[eeg_data] = ft_preprocessing(cfg,eeg_data);
%% average EOGs
cfg = [];
cfg.avgoverchan ='yes';
[eeg_data] = ft_selectdata(cfg,eeg_data);
%% concatenate EOG-data across trials
nsamp = zeros(length(eeg_data.trial),1);
for kt = 1:length(eeg_data.trial)
    nsamp(kt) = length(eeg_data.trial{kt});
end;
concat = zeros(length(eeg_data.label),sum(nsamp));
for jt = 1:length(eeg_data.label)
    idx = 1:length(eeg_data.trial{1});
    for kt = 1:length(eeg_data.trial)
        concat(jt,idx) = eeg_data.trial{kt}(jt,:);
        if sign(kt-length(eeg_data.trial))==-1
            idx = idx(end)+1:idx(end)+length(eeg_data.trial{kt+1});
        end;
    end;
    concat(jt,:) = concat(jt,:) - mean(concat(jt,:),2);
    concat(jt,:) = concat(jt,:)./std(concat(jt,:),0,2);
end;
%% concatenate MEG-data across trials
nsamp = zeros(length(comp_data.trial),1);
for kt = 1:length(comp_data.trial)
    nsamp(kt) = length(comp_data.trial{kt});
end;
concat2 = zeros(length(comp_data.label),sum(nsamp));
for jt = 1:length(comp_data.label)
    idx = 1:length(comp_data.trial{1});
    for kt = 1:length(comp_data.trial)
        concat2(jt,idx) = comp_data.trial{kt}(jt,:);
        if kt <length(comp_data.trial)
            idx = idx(end)+1:idx(end)+length(comp_data.trial{kt+1});
        end;
    end;
    concat2(jt,:) = concat2(jt,:) - mean(concat2(jt,:),2);
    concat2(jt,:) = concat2(jt,:)./std(concat2(jt,:),0,2);
end;
%% search EOG data for significant events (treshold based)
trsh = 2;
idx = cell(size(concat,1),1);
for jt = 1:size(concat,1)
    if sum(concat(jt,:) >trsh,2) > sum(concat(jt,:) <-trsh,2)
        x = concat(jt,:).*(concat(jt,:)>trsh);
        x = x >median(x(x~=0));
    else
        x = concat(jt,:).*(concat(jt,:)<-trsh);
        x = x < median(x(x~=0));
    end;
    idx{jt} = find(sign(diff(x))==1)+1;
end;
%%
if ~isempty(idx)
    %% align both EOG and MEG data around significant EOG events
    protoBlink = cell(length(idx),1);
    bcComp = cell(length(idx),1);
    k = 0;
    for jt = 1:length(idx)
        protoBlink{jt} = zeros(length(idx{jt}),comp_data.fsample*1+1);
        bcComp{jt} = zeros(length(comp_data.label),length(idx{jt}),comp_data.fsample*1+1);
        for kt = 1:length(idx{jt})
            if (idx{jt}(kt)-0.5*comp_data.fsample>0) && (idx{jt}(kt)+0.5*comp_data.fsample<size(concat,2))
                k = k+1;
                protoBlink{jt}(kt,:) = concat(jt,idx{jt}(kt)-0.5*comp_data.fsample:idx{jt}(kt)+0.5*comp_data.fsample);
                bcComp{jt}(:,kt,:) = concat2(:,idx{jt}(kt)-0.5*comp_data.fsample:idx{jt}(kt)+0.5*comp_data.fsample);
            end;
        end;
    end;
    %%
    if ~isempty(idx{1})
        dum1 = struct;
        dum1.label(1) = {'protoBink1'};
        dum1.trial = cell(1,size(protoBlink{1},1));
        dum1.time = cell(1,size(protoBlink{1},1));
        for nt = 1:size(protoBlink{1},1)
            dum1.trial{nt} = protoBlink{1}(nt,:);
            dum1.time{nt} = 0:1/eeg_data.fsample:(size(protoBlink{1},2)-1)/eeg_data.fsample;
        end;
        dum1 = ft_timelockanalysis([],dum1);
        protoBlink{1} = dum1.avg;
        
    end;
    if ~isempty(idx{1})
        dum3 = struct;
        dum3.label = comp_data.label;
        dum3.trial = cell(1,size(bcComp{1},2));
        dum3.time = cell(1,size(bcComp{1},2));
        for nt = 1:size(bcComp{1},2)
            dum3.trial{nt} = squeeze(bcComp{1}(:,nt,:));
            dum3.time{nt} = 0:1/eeg_data.fsample:(size(bcComp{1},3)-1)/eeg_data.fsample;
        end;
        dum3 = ft_timelockanalysis([],dum3);
        bcComp{1} = dum3.avg;
        
    end;
    
    try
        if ~isempty(idx{2})
            dum2 = struct;
            dum2.label(1) = {'protoBink2'};
            dum2.trial = cell(1,size(protoBlink{2},1));
            dum2.time = cell(1,size(protoBlink{2},1));
            for nt = 1:size(protoBlink{2},1)
                dum2.trial{nt} = protoBlink{2}(nt,:);
                dum2.time{nt} = 0:1/eeg_data.fsample:(size(protoBlink{2},2)-1)/eeg_data.fsample;
            end;
            dum2 = ft_timelockanalysis([],dum2);
            protoBlink{2} = dum2.avg;
            
        end;
        
        if ~isempty(idx{2})
            dum4 = struct;
            dum4.label = comp_data.label;
            dum4.trial = cell(1,size(bcComp{2},2));
            dum4.time = cell(1,size(bcComp{2},2));
            for nt = 1:size(bcComp{2},2)
                dum4.trial{nt} = squeeze(bcComp{2}(:,nt,:));
                dum4.time{nt} = 0:1/eeg_data.fsample:(size(bcComp{2},3)-1)/eeg_data.fsample;
            end;
            dum4 = ft_timelockanalysis([],dum4);
            bcComp{2} = dum4.avg;
            
        end;
    catch
        %do nothing;
    end;   
    %%
    idx2 = zeros(length(idx),1);
    for kt = 1:length(idx)
        idx2(kt) = ~isempty(idx{kt});
    end;
    idx3 = find(idx2);
    %%
    if ~isempty(idx3)
        rxy = cell(length(idx3),1);
        for kt = 1:length(idx3)
            rxy{kt}.r = zeros(length(comp_data.label),1);
            
            for jt = 1:length(comp_data.label);
                rxy{kt}.r(jt) = max(abs(xcorr(protoBlink{idx3(kt)},bcComp{idx3(kt)}(jt,:),'coeff')));
            end;
        end;
        
        idx = [];
        for kt = 1:length(rxy)
            idx = [idx;find(rxy{kt}.r >=.7)];
        end;
        SSpidx = unique(idx);
        
    end;  
    
    tx= -.5:1/comp_data.fsample:.5;
    t_idx = find(tx==0);
    
    m = cell(length(bcComp),1);
    for kt = 1:length(bcComp)
        for mt = 1:length(SSpidx)
            M = mean(bcComp{kt}(SSpidx(mt),:),2);
            SD = std(bcComp{kt}(SSpidx(mt),:),[],2);
            bcComp{kt}(SSpidx(mt),:) = (bcComp{kt}(SSpidx(mt),:) - M)./SD;
            m{kt}(mt) = max(bcComp{kt}(SSpidx(mt),t_idx));
        end;
    end;
    
    sel_idx = [];
    for kt = 1:length(m)
        sel_idx = [sel_idx;find(sign(m{kt}-2)==1)];
    end;
    sel_idx = unique(sel_idx);
    
    SSpidx = SSpidx(sel_idx);
end;