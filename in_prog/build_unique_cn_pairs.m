function [stim_mat] = build_unique_cn_pairs(ID,params)
flag = 0;
while flag <1
    %%
    stim_mat.ID = ID.id;
    stim_mat.idx =ID.idx;
    %%
    ID2 = ID.id;
    %%
    %ID2 = cell(length(ID),1);
    idx = zeros(length(ID2),1);
    for it = 1:length(ID2)
        
        %ID2{it} = ID(it)(1:regexp(ID(it),'.jpg')-1);
        idx(it) = it;
        
    end;
    %%
    concept_n = params.concept_neurons;
    cc = zeros(length(idx),1);
    for it = 1:length(concept_n)
        
        chck = regexp(ID2,concept_n{it});
        
        chck2 = zeros(length(chck),1);
        for jt = 1:length(chck)
            chck2(jt) = isempty(chck{jt});
        end;
        
        ix = find(chck2 ==0);
        cc(ix) = it+zeros(length(ix),1);
        
    end;
    %%
    n = length(concept_n);
    
    cp = zeros(n*(n-1)/2+n,2);
    k = 0;
    
    x1 = 1;
    for it = 1:n
        for jt = x1:n
            k=k+1;
            cp(k,:) = [it jt];
            
        end;
        x1 = x1+1;
    end;
    %%
    samp = [idx cc];
    
    [c_idx] = find(samp(:,2)==0);
    [c_idx] = c_idx(randperm(length(c_idx)));
    
    trl = 0;
    seq = zeros(3,round(length(find(samp(:,2)~=0))/2));
    f1 = 0;
    h = [];
    while f1==0
        for it = 1:size(cp,1)
            
            p = cp(it,:);
            
            sel = [];
            ix = cell(length(p),1);
            for jt = 1:length(p)
                ix{jt} = find(samp(:,2) == p(jt));
            end;
            if isempty(ix{1})
                f1=1;
            else
                ix{1} = ix{1}(randperm(length(ix{1})));
                sel(1) = ix{1}(1);
                
                
                if (length(find(ismember(samp(:,2),cp(it,:)))) >=2) && ~isempty(ix{2})
                    f2 = 0;
                    while f2 ==0
                        ix{2} = ix{2}(randperm(length(ix{2})));
                        sel(2) = ix{2}(1);
                        if sel(2) ~= sel(1)
                            f2=1;
                        end;
                    end;
                end;
                
                if length(sel) ==2
                    trl = trl+1;
                    h(trl,:) = p;
                    seq(2:3,trl) = samp(sel,1);
                    samp(sel,:) = [];
                end;
            end;
        end;
        
    end;
    %%
    %%
    for it = 1:size(seq,2)
        seq(1,it) = c_idx(it);
    end;
    %%
    [~,d2] = find(seq(2:3,:) == 0);
    
    rand('state',sum(100*clock));
    d2 = d2(randperm(length(d2)));
    d2 = d2(randperm(length(d2)));
    d2 = d2(randperm(length(d2)));
    
    seq(:,d2) = [];
    %% sanity check
    x = seq(2:3,:);
    x = x(:);
    if (length(x) == length(unique(x))) ~= 1
        error('non-unique event codes not permitted');
    end;
    %%
    stim_mat.tc = [];
    stim_mat.lkp = [];
    stim_mat.c = find(cc==0);
    stim_mat.p = find(cc ~=0);
    stim_mat.seq = seq;
    stim_mat.xc =[];
    if any(ismember(stim_mat.c,stim_mat.p))
        error('overlap between pair and cue indices');
    end;
    
    c =0;
    chck = [];
    for kt = 1:size(stim_mat.seq,2)
        sel = stim_mat.seq(2:3,kt);
        for lt = 1:size(params.stim_mat.seq,2)
            sel2 = params.stim_mat.seq(2:3,lt);
            if isequal(sel,sel2)
                c =c+1;
                chck(c) = 1;
                fprintf('warning searching for compatible cn-pair configuration\n');
            end;
        end;
    end;
    
    if sum(chck)==0
        flag=1;
    end;
end;
fprintf('exiting cn pair configuration: all pairs ok\n');
return;