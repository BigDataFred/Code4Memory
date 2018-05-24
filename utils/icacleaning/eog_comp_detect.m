function [eog_indx] = eog_comp_detect(eog_data,comp,template)

%% detection of ICs with Blink artifacts
% compute correlation of time courses between ICs and vertical EOG signal
%%
ts= zeros(length(comp.time),1);for it = 1:length(comp.time);ts(it) = length(comp.time{it});end;ts = sum(ts); 
concat = zeros(length(comp.topolabel),ts);
indx = 1:length(comp.time{1});
for it = 1:length(comp.trial)
    
    
    concat(:,indx) =comp.trial{it};
    if it < length(comp.trial)
        indx = indx(end)+1:indx(end)+length(comp.time{it+1});
    end;
end;
%%
ts= zeros(length(eog_data.time),1);for it = 1:length(eog_data.time);ts(it) = length(eog_data.time{it});end;ts = sum(ts); 
eog_concat = zeros(length(eog_data.label),ts);
indx = 1:length(eog_data.time{1});
for it = 1:length(eog_data.trial)
        
    eog_concat(1,indx) = eog_data.trial{it}(1,:);
    eog_concat(2,indx) = eog_data.trial{it}(2,:);
    
    if it < length(eog_data.trial)
        indx = indx(end)+1:indx(end)+length(eog_data.time{it+1});
    end;
    
end;

%%
dummy = struct;
dummy.fsample = comp.fsample;
dummy.label = comp.topolabel;
dummy.trial{1} = concat;
dummy.time{1} = [1:size(concat,2)]./comp.fsample;
dummy.cfg = [];

cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = [1 15];
cfg.bpfilttype = 'but';
cfg.bpfiltdir = 'twopass';
cfg.bpfiltorder = 4;%6*fix(dummy.fsample/cfg.bpfreq(1));

filt_concat = ft_preprocessing(cfg,dummy);

dummy.label = eog_data.label;
dummy.trial{1} = eog_concat;
filt_eog = ft_preprocessing(cfg,dummy);


r1 = zeros(length(comp.topolabel),1);
r2 = zeros(length(comp.topolabel),1);
p1 = zeros(length(comp.topolabel),1);
p2 = zeros(length(comp.topolabel),1);
for jt = 1:length(comp.topolabel)
    [r1(jt),p1(jt)] = corr(filt_eog.trial{1}(1,:)',filt_concat.trial{1}(jt,:)');% corr with vEOG
    [r2(jt),p2(jt)] = corr(filt_eog.trial{1}(2,:)',filt_concat.trial{1}(jt,:)');% corr with hEOG
end;
indx1 = find(r1>=.5);
indx2 = find(r2>=.5);

if size(indx1,2) > size(indx1,1)
    indx1 = indx1';
    indx2 = indx2';
end;
eog_indx1 = unique([indx1;indx2]);
%% match IC topography with template EOG topography
[temp_topo] = eog_topography_template(template);
%%
r1 = zeros(length(comp.topolabel),1);
r2 = zeros(length(comp.topolabel),1);
r3 = zeros(length(comp.topolabel),1);
r4 = zeros(length(comp.topolabel),1);
r5 = zeros(length(comp.topolabel),1);
r6 = zeros(length(comp.topolabel),1);
r7 = zeros(length(comp.topolabel),1);
for jt = 1:length(comp.topolabel)
    [r1(jt),p1] = corr(temp_topo.topo1,comp.topo(:,jt),'Type','Spearman');
    [r2(jt),p2] = corr(temp_topo.topo2,comp.topo(:,jt),'Type','Spearman');
    [r3(jt),p3] = corr(temp_topo.topo3,comp.topo(:,jt),'Type','Spearman');
    [r4(jt),p4] = corr(temp_topo.topo4,comp.topo(:,jt),'Type','Spearman');
    [r5(jt),p5] = corr(temp_topo.topo5,comp.topo(:,jt),'Type','Spearman');
    [r6(jt),p6] = corr(temp_topo.topo6,comp.topo(:,jt),'Type','Spearman');
    [r7(jt),p7] = corr(temp_topo.topo7,comp.topo(:,jt),'Type','Spearman');
end;
eog_indx2 =unique([find(r1>=0.5);find(r2>=0.5);find(r3<=0.5);find(r4<=0.5);find(r5<=0.5);find(r6<=0.5);find(r7<=0.5)]);
if size(eog_indx1,2) > size(eog_indx1,1)
    eog_indx1 = eog_indx1';
    eog_indx1 = eog_indx1';
end;
if size(eog_indx2,2) > size(eog_indx2,1)
    eog_indx2 = eog_indx2';
    eog_indx2 = eog_indx2';
end;
eog_indx = intersect(eog_indx1,eog_indx2);