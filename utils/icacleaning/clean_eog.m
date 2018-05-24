function [eog_idx,eog,SSp] = clean_eog(ica_data,eog_data)
%%
addpath('~froux/froux/sl_meg/mcode/icacleaning/');
%%
[eog_idx] = eog_comp_detect2(ica_data,eog_data);
%%
if ~isempty(eog_idx)
    cfg = [];
    cfg.eog_idx = eog_idx;
    cfg.bw = 0.15;
    
    [eog_idx2,~] = eog_comp_detect3(cfg,ica_data);
    
    if size(eog_idx2,2)>size(eog_idx2,1)
        eog_idx2 = eog_idx2';
    end;
    eog_idx = sort([eog_idx;eog_idx2]);        

end;
%%
[SSpidx,protoBlink,bcComp] = SSp_comp_detect(eog_data,ica_data);

SSp.idx = SSpidx;
SSp.protoBlink = protoBlink;
SSp.bcComp = bcComp;
%%
if size(SSpidx,2)>size(SSpidx,1)
    SSpidx = SSpidx';
end;
%%
eog_idx = unique(sort([eog_idx;SSpidx]));

cfg =[];
cfg.eog_idx = eog_idx;
cfg.bw = 0.1;

[~,bData] = eog_comp_detect3(cfg,ica_data);

eog.idx = eog_idx;
eog.bData = bData;
