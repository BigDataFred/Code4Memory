%%
addpath(genpath('/media/rouxf/rds-share/Fred/code/mcode/custom/'));
addpath(genpath('/media/rouxf/rds-share/Fred/code/mcode/toolboxes/osort-v3-rel/'));

%%
rpath = '/media/rouxf/rds-share/iEEG_DATA/MICRO/P02/fVSpEM/';
load([rpath,'Cluster_info_P02.mat']);

chck = dir(rpath);
chck(1:2) = [];

sesh_id = {chck(find([chck(:).isdir])).name};

%%
for it = 1%:length(sesh_id)
    
    %search for files with raw spike data
    CSCfiles = dir([rpath,sesh_id{it},filesep,'spike_dat',filesep,'*.mat']);
    
    % extract index for current session
    sesh_ix = find(clusterInfo(:,1)==it);
    
    % extract index for microwires with putative SUs
    mw_ix = clusterInfo(sesh_ix,2);
    
    % only keep MWs with SUs
    CSCfiles = CSCfiles(unique(mw_ix));
    
    % IDs of MWs
    mwID = unique(clusterInfo(sesh_ix,find(strcmp(readme,'microwireID'))));
    
    if length(CSCfiles) ~= length(mwID) %n1 vs n2
        error('n1 and n2 must have equal length');
    end;
    
    count = 0;
    SU = cell(1,length(mwID));
    spk_data2 = {};

    for jt = 1:length(mwID)% loop over MWs
        
        %load data for each MW
        dat = load([rpath,sesh_id{it},filesep,'spike_dat',filesep,CSCfiles(jt).name]);
        
        spk_dat = dat.save_data{1}{1}{1};        
                        
        %%

        %indexes of clusters on each MW
        clu_ix = clusterInfo(sesh_ix(find(clusterInfo(sesh_ix,2) == mwID(jt))),3);
        
        %%
        dum = [];
        ix = [];
        for kt = 1:length(spk_dat.waveform)
            dum = [dum;squeeze(spk_dat.waveform{kt})'];
            ix(kt) = size(squeeze(spk_dat.waveform{kt})',1);
        end;
        [~,score] = pca( dum );        
        pcs = [score(:,1) score(:,2) score(:,3)];

        makePrinComp_plot(pcs,ix,spk_dat.waveform,clu_ix)
        
%         %%
%         spike_stats = {};
%         for kt = 1:length( spk_dat.waveform )
%             timestamps = spk_dat.timestamp{kt};
%             spike_stats{kt} = getStatsForCluster([],timestamps);            
%         end;
%         plotClusterStats(spike_stats);
%         
        %%
        txt = [];
        for kt = 1:length(clu_ix)
            txt = [txt [num2str(clu_ix(kt)),',']];
        end;
        txt(end) = [];
        
        [resp] = input(['Enter number of bad cluster(s): options are ',txt],'s');
        for kt = 1:length(resp)
            resp(kt) = str2double(resp(kt));
        end;
        del_ix = find(ismember(clu_ix,resp));
        del_ix2 = find(clusterInfo(sesh_ix,2) == mwID(jt));
        del_ix2 = del_ix2(del_ix);
        
        clusterInfo(sesh_ix(del_ix2),:) = [];
        sesh_ix(del_ix2) = [];
        clu_ix(del_ix) = [];
        
        if ~isempty(resp)
            dum = [];
            ix = [];
            for kt = 1:length(spk_dat.waveform)
                dum = [dum;squeeze(spk_dat.waveform{kt})'];
                ix(kt) = size(squeeze(spk_dat.waveform{kt})',1);
            end;
            [~,score] = pca( dum );
            pcs = [score(:,1) score(:,2) score(:,3)];
            
            makePrinComp_plot(pcs,ix,spk_dat.waveform,clu_ix)
            
        end;
        close all;
        
        %%
        count = count+1;
        spk_dat2{count} =  extract_cluster_field(spk_dat,clu_ix);
        
        % if there are more than 2 cluster do projection test
        if length( clu_ix ) > 1
            
            p = [];
            c=0;
            for kt = 1:length(clu_ix)
                for lt = kt+1:length(clu_ix)
                    c = c+1;
                    p(c,:) = [kt lt];
                end;
            end;
            
            ix = zeros(1,size(p,1));
            for kt = size(p,1)%1:size(p,1)
                
                spikes = cell(1,size(p,2));
                for lt = 1:size(p,2)
                    spikes{lt} = squeeze( spk_dat2{count}.waveform{p(kt,lt)} )';                    
                end;
                
                [~,~, ~,~,~,d] = projectionTest( spikes{1},spikes{2} );
                
                if d < 5
                    ix(kt) = 1;
                    return;
                end;
                
            end;
            
            merge_ix = find(ix ==1);
            
            dum = [];
            for kt = 1:length(merge_ix)
                return;
                ix = p(merge_ix(kt),:);
                if isempty(dum)
                    dum = merge_clusters([],spk_dat2,ix);
                else
                    dum = merge_clusters(dum,spk_dat2,ix);
                end;
                spk_dat2 =  remove_cluster_field(spk_dat2,ix(1));
                spk_dat2 =  remove_cluster_field(spk_dat2,ix(2));
            end;                           
        end;
        
%         nCluster = length(getfield(spk_dat,'label'));
%         for kt = 1:nCluster
%             count = count+1;
%             [SU{count}] = extract_SU_parameters(spk_dat,kt);
%         end;
        
    end;
    
end;