function [sortedSpikes,wvltCoeffs] = doSpikeSorting_waveclus( waveclusdata , par )

sortedSpikes = [];
[num_temp] = floor((par.maxtemp-par.mintemp)/par.tempstep); 

naux = min(par.max_spk,size(waveclusdata.spikes,1));
par.min_clus = max(par.min_clus_abs,par.min_clus_rel*naux);

inspk = wave_features(waveclusdata.spikes , par );

if strcmp(par.permut,'n')
    
    if size(waveclusdata.spikes,1) > par.max_spk % goes for template matching if too many spikes
        inspk_aux = inspk(1:naux,:);
    else
        inspk_aux = inspk;
    end;
    
    save(par.fname_in,'inspk_aux','-ascii');
    [clu,tree] = run_cluster( par );
    
    size(tree)
    if ~isempty( tree ) && ( size(tree,1)-1 >= num_temp-1 )
        [temp] =  find_temp(tree , par);
        
        %define clusters
        c = {};
        c{1} = find(clu(temp,3:end)==0);
        c{2} = find(clu(temp,3:end)==1);
        c{3} = find(clu(temp,3:end)==2);
        c{4} = find(clu(temp,3:end)==3);
        c{5} = find(clu(temp,3:end)==4);
        c{6} = setdiff(1:size(waveclusdata.spikes,1), sort([c{1} c{2} c{3} c{4} c{5}]));
    end;
else
    
    if size(waveclusdata.spikes,1) > par.max_spk % goes for template matching if too many spikes
        iperm = randperm(length(inspk));
        iperm(naux+1:end) = [];
        inspk_aux = inspk(iperm,:);
    else
        iperm = randperm(length(inspk));
        inspk_aux = inspk(iperm,:);
    end;
    
    save(par.fname_in,'inspk_aux','-ascii');
        
    [clu,tree] = run_cluster( par );
    size(tree)
    if ~isempty(tree) && ( size(tree,1)-1 >= num_temp-1 )

        %try
        [temp] =  find_temp(tree , par);
        %catch
        %return
        %end;
        
        c = {};
        c{1} = iperm(find(clu(temp,3:end)==0));
        c{2} = iperm(find(clu(temp,3:end)==1));
        c{3} = iperm(find(clu(temp,3:end)==2));
        c{4} = iperm(find(clu(temp,3:end)==3));
        c{5} = iperm(find(clu(temp,3:end)==4));
        c{6} = setdiff(1:size(waveclusdata.spikes,1), sort([c{1} c{2} c{3} c{4} c{5}]));
    end;
end;

% it template matching was done then force
if (size(waveclusdata.spikes,1) > par.max_spk) || (strcmp(par.foce_auto,'auto'))
    
    cl = zeros(size(waveclusdata.spikes,1),1);
    if length(c{1}) >=par.min_clus;cl(c{1}) = 1; end;
    if length(c{2}) >=par.min_clus;cl(c{2}) = 2; end;
    if length(c{3}) >=par.min_clus;cl(c{3}) = 3; end;
    if length(c{4}) >=par.min_clus;cl(c{4}) = 4; end;
    if length(c{5}) >=par.min_clus;cl(c{5}) = 5; end;
    
    f_in = waveclusdata.spikes(cl ~=0,:);
    f_out = waveclusdata.spikes(cl ==0,:);
    cl_in = cl(find(cl ~=0),:);
    cl_out = force_membership_wc(f_in, cl_in, f_out, par);
    cl(cl==0) = cl_out;
    c = {};
    c{6} = find(cl==0);
    c{1} = find(cl==1);
    c{2} = find(cl==2);
    c{3} = find(cl==3);
    c{4} = find(cl==4);
    c{5} = find(cl==5);
    
end;

cluster = zeros(size(waveclusdata.spikes,1),2);
cluster(:,2) = waveclusdata.index;

num_cluster = length(find([length(c{1}) length(c{2}) length(c{3}) length(c{4}) length(c{5}) length(c{6})]) >= par.min_clus);

clus_pop = [];
clus_pop = [clus_pop length([c{6}])];

for it = 1:length(c)-1
    if length(c{it}) > par.min_clus
        
        clus_pop = [clus_pop length([c{it}])];
        cluster(c{it},1) = it;
        
    end;
end;

sortedSpikes.newSpikeTimes = cluster(:,2)';
sortedSpikes.assignedCluster = cluster(:,1)';
sortedSpikes.wavf = waveclusdata.spikes;
sortedSpikes.num_clus = length(clus_pop)-1;

wvltCoeffs = inspk;

return;
