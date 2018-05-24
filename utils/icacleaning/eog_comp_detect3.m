function [eog_idx,bData] = eog_comp_detect3(cfg,comp_data)
if isfield(cfg,'eog_idx')
    cfg.eog_idx = unique(cfg.eog_idx);
    %%
    concat = zeros(length(comp_data.label),length(comp_data.trial)*length(comp_data.trial{1}));
    z = zeros(length(comp_data.label),length(comp_data.trial)*length(comp_data.trial{1}));

    for jt = 1:length(comp_data.label)
        idx = 1:length(comp_data.trial{1});
        for kt = 1:length(comp_data.trial)
            concat(jt,idx) = comp_data.trial{kt}(jt,:);
            if kt < length(comp_data.trial)
                idx = idx(end)+1:idx(end)+length(comp_data.trial{kt+1});
            end;
        end;
        z(jt,:) = concat(jt,:) - mean(concat(jt,:),2);
        z(jt,:) = concat(jt,:)./std(concat(jt,:),0,2);
    end;
    %%
    idx = cell(1,length(cfg.eog_idx));
    trsh = 2;
    for it = 1:length(cfg.eog_idx)
%         dum             = struct;
%         dum             = struct;
%         dum.fsample     = comp_data.fsample;
%         dum.label       = comp_data.label;
%         dum.trial{1}    = concat;
%         dum.time{1}     = 0:1/dum.fsample:(length(dum.trial{1})-1)/dum.fsample;
%         dum.cfg         = [];
%         
%         cfg2            = [];
%         cfg2.trl        = [comp_data.sampleinfo zeros(size(comp_data.sampleinfo,1),1)];
%         
%         dum = ft_redefinetrial(cfg2,dum);
%         
%         cfg2                            = [];
%         cfg2.trl                        = [dum.sampleinfo zeros(size(dum.sampleinfo,1),1)];
%         cfg2.continuous                 = 'yes';
%         cfg2.artfctdef.eog.cutoff       =4;
%         cfg2.artfctdef.eog.channel      = dum.label{cfg.eog_idx};
%         cfg2.artfctdef.eog.interactive =  'no';
%         cfg2.artfctdef.eog.inspect     = dum.label{cfg.eog_idx};
%         cfg2.artfctdef.eog.trlpadding   = 0;
%         cfg2.artfctdef.eog.fltpadding   = 0;
%         cfg2.artfctdef.eog.artpadding   = 0;
%         
%         [cfg2,eog_art]  = ft_artifact_eog(cfg2,dum);
        
        [eog_sig] = concat(cfg.eog_idx(it),:);
        eog_sig = conv(eog_sig,gausswin(100),'same');
        z = (eog_sig - mean(eog_sig))./std(eog_sig);
        
        [~,idx{it}] = localmax(z);
        idx{it} = idx{it}(find(sign(z(idx{it})-trsh)==1));
        idx{it}(find(sign(diff(idx{it})-comp_data.fsample*cfg.bw)==-1)+1) = [];
    end;
    %%
%     trsh = 2;
%     x = zeros(size(z));
%     idx = cell(1,length(cfg.eog_idx));
%     for it = 1:length(cfg.eog_idx)
%         if sum(z(cfg.eog_idx(it),:) >trsh,2) > sum(z(cfg.eog_idx(it),:) <-trsh,2)
%             x(it,:) = z(cfg.eog_idx(it),:).*(z(cfg.eog_idx(it),:)>trsh);
%             x(it,:) = x(it,:) >max(x(it,(x(it,:)~=0)))/4;
%         else
%             x(it,:) = z(cfg.eog_idx(it),:).*(z(cfg.eog_idx(it),:)<-trsh);
%             x(it,:) = x(it,:) < min(x(it,(x(it,:)~=0)))/4;
%         end;
%         idx{it} = find(sign(diff(x(it,:)))==1)+1;
%     end;

    % %%
    % figure;
    % subplot(211);
    % a = gca;
    % plot(concat(cfg.eog_idx,:));
    % subplot(212);
    % a(2) = gca;
    % hold on;
    % plot(1:length(x),x);
    % plot(idx,x(idx),'r*');
    % ylim([-2 2]);
    %
    % set(a,'Xlim',[1*comp_data.fsample 125*comp_data.fsample]);
    %%
    protoBlink  = cell(length(cfg.eog_idx),1);
    bcMEG       = cell(length(cfg.eog_idx),1);
    avg1        = cell(length(cfg.eog_idx),1);
    avg2        = cell(length(cfg.eog_idx),1);
    zavg1        = cell(length(cfg.eog_idx),1);
    zavg2        = cell(length(cfg.eog_idx),1);
    
    idx2 = [];
    idx3 = [];
    %%
    for it = 1:length(cfg.eog_idx)
        protoBlink{it}.trl = zeros(length(idx{it}),comp_data.fsample*cfg.bw*2+1);
        bcMEG{it}.trl = zeros(length(comp_data.label),length(idx{it}),comp_data.fsample*cfg.bw*2+1);
        k = 0;
        protoBlink{it}.sampleinfo = zeros(length(idx{it}),2);
        for jt = 1:length(idx{it})
            if (idx{it}(jt)-(cfg.bw*comp_data.fsample)>0) && (idx{it}(jt)+(cfg.bw*comp_data.fsample)<size(concat,2))
                k = k+1;
                protoBlink{it}.trl(k,:) = concat(cfg.eog_idx(it),idx{it}(jt)-(cfg.bw*comp_data.fsample):idx{it}(jt)+(cfg.bw*comp_data.fsample));
                bcMEG{it}.trl(:,k,:) = concat(:,idx{it}(jt)-(cfg.bw*comp_data.fsample):idx{it}(jt)+(cfg.bw*comp_data.fsample));
                
                protoBlink{it}.sampleinfo(jt,:) = [idx{it}(jt)-(cfg.bw*comp_data.fsample) idx{it}(jt)+(cfg.bw*comp_data.fsample)];
                
            end;
        end;
        %%
        avg1{it} = mean(protoBlink{it}.trl,1);
        avg2{it} = squeeze(mean(bcMEG{it}.trl,2));
        r = zeros(size(avg2,1),1);
        for kt = 1:size(avg2)
            
            r(kt) = max(abs(xcorr(avg2{it}(kt,:),avg1{it},'coeff')));
            
        end;
        
        stp = find(sign(r-.9)==1);
        if ~isempty(stp)
            if sign(size(stp,1)-size(stp,2))==1
                stp = stp';
            end;
            idx2 = [idx2 stp];
        end;
        
        M = repmat(mean(avg1{it},2),[1 size(avg1{it},2)]);
        SD = repmat(std(avg1{it},[],2),[1 size(avg1{it},2)]);
        
        zavg1{it} = avg1{it} -M;
        zavg1{it} = zavg1{it}./SD;
        
        M = repmat(mean(avg2{it},2),[1 size(avg2{it},2)]);
        SD = repmat(std(avg2{it},[],2),[1 size(avg2{it},2)]);
        
        zavg2{it} = avg2{it} -M;
        zavg2{it} = zavg2{it}./SD;
        
        if size(idx3,2) > size(idx3,1)
            idx3 = idx3';
        end;
        
        idx3 = [idx3;find(sign(max(abs(zavg2{it}),[],2)-2)==1)];
        
    end;
    %%
    bData = cell(length(protoBlink),1);
    for it =1:length(protoBlink)
        bData{it}.trl = protoBlink{it}.trl;
        bData{it}.MEGtrl = bcMEG{it};
        bData{it}.sampleinfo = protoBlink{it}.sampleinfo;
        bData{it}.avg = avg1{it};
        bData{it}.MEGavg = avg2{it};
        bData{it}.zavg = zavg1{it};
        bData{it}.zMEGavg = zavg2{it};
    end;
    %%
    [eog_idx] = intersect(idx2,idx3);
    % %%
    % cfg = [];
    % cfg.layout = 'CTF275.lay';
    % cfg.parameter = 'avg';
    %
    % figure;
    % dum = struct;
    % dum.time = 0;
    % dum.avg = comp_data.topo(:,cfg.eog_idx);
    % dum.dimord = 'chan_time';
    % dum.label = meg_data.label;
    %
    % ft_topoplotER(cfg,dum);
    %
    %
    % figure;
    % subplot(5,1,1:4);
    % imagesc(linspace(-.5,.5,201),1:length(idx),abs(protoBlink));
    % axis xy;
    % subplot(515);
    % hold on;
    % plot(linspace(-.5,.5,201),abs(mean(protoBlink,1)));
    % %%
    % cfg = [];
    % cfg.layout = 'CTF275.lay';
    % cfg.parameter = 'avg';
    %
    % figure;
    % dum = struct;
    % dum.time = 0;
    % dum.avg = comp_data.topo(:,idx);
    % dum.dimord = 'chan_time';
    % dum.label = meg_data.label;
    %
    % ft_topoplotER(cfg,dum);
    %
    %
    % figure;
    % subplot(5,1,1:4);
    % imagesc(linspace(-.5,.5,201),1:length(idx),(abs(squeeze(bcMEG(idx(1),:,:)))));
    % axis xy;
    % subplot(515);
    % hold on;
    % plot(linspace(-.5,.5,201),abs(bcMEG2(idx(1),:)));
    %%
    %cfg = [];cfg.layout = 'CTF275.lay';cfg.viewmode = 'comp_dataonent';cfg.channel = [idx cfg.eog_idx];ft_databrowser(cfg,comp_data);
else
    eog_idx = [];
end;