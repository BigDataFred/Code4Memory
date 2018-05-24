%% include in path
addpath(genpath('/media/rouxf/rds-share/Fred/code/mcode/custom/'));
addpath(genpath('/media/rouxf/rds-share/Fred/code/mcode/toolboxes/wave_clus-testing/'));

%% IDs of patients to include
pID = {'P02' 'P04' 'P05' 'P22AMS' 'P23AMS'};%patient ID

%% controls the binning of the average spike-density
step = 500;

%%
FRce = [];
FRc = [];
FRe = [];
nItemsCE = [];
nItemsC = [];
nItemsE = [];
M1 = [];
M2 = [];
M3 = [];

c1 = 0;
c2 = 0;
c3 = 0;

%% loop over patients
for pt = 1:length(pID)
    
    fprintf([num2str(pt),'/',num2str( length(pID) ),'\n']);
    
    %% set the path for data access
    rpath = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',pID{pt},'/'];%rooth path
    
    [exp,~] = get_EM_rec_info(rpath);% extract session labels
    
    [rpath] = [rpath,exp,filesep]; % path to read data
    
    %%
    pooledDat = dir([rpath,pID{pt},'_pooledSPKdataEMtask.mat']);
    load([rpath,pooledDat.name]);
    
    %%
    dt2 = dat.dt(1):step:dat.dt(end);
    for ct = 1:length(dat.cluRaster)
        
        selIdx1 = dat.selIdxCE{ct};
        selIdx2 = dat.selIdxC{ct};
        selIdx3 = dat.selIdxE{ct};
        
        for gt = 1:length(selIdx1)
            
            n = [];
            for it = 1:length( dt2 )-1
                idx = find(dat.dt >=dt2(it) & dat.dt < dt2(it+1));
                n(:,it) = sum(dat.cluRaster{ct}{selIdx1(gt)}(:,idx),2);
            end;
            
            n = n./(step/dat.Fs);
            %n(:,[1 end]) = [];
            
            M = mean(mean(n(:,find(dt2 <0)),2));
            SD = std(mean(n(:,find(dt2 <0)),2),0,1);
            M = repmat(M,size(n));
            SD = repmat(SD,size(n));
            n = (n-M)./SD;
            n(isnan(n)) = 0;
            n(isinf(n)) = 0;
            
            if any(isnan(n))
                error('NaNs detected');
            end;
            
            c1=c1+1;
            M1(c1,:) = sum(n,1)./size(n,1);
            
        end;
        
        for gt = 1:length(selIdx2)
            
            n = [];
            for it = 1:length( dt2 )-1
                idx = find(dat.dt >=dt2(it) & dat.dt < dt2(it+1));
                n(:,it) = sum(dat.cluRaster{ct}{selIdx2(gt)}(:,idx),2);
            end;
            
            n = n./(step/dat.Fs);
            %n(:,[1 end]) = [];
            
            M = mean(mean(n(:,find(dt2 <0)),2));
            SD = std(mean(n(:,find(dt2 <0)),2),0,1);
            M = repmat(M,size(n));
            SD = repmat(SD,size(n));
            n = (n-M)./SD;
            n(isnan(n)) = 0;
            n(isinf(n)) = 0;
            
            if any(isnan(n))
                error('NaNs detected');
            end;
            
            c2=c2+1;
            M2(c2,:) = sum(n,1)./size(n,1);
            
        end;
        
        for gt = 1:length(selIdx3)
            
            n = [];
            for it = 1:length( dt2 )-1
                idx = find(dat.dt >=dt2(it) & dat.dt < dt2(it+1));
                n(:,it) = sum(dat.cluRaster{ct}{selIdx3(gt)}(:,idx),2);
            end;
            
            n = n./(step/dat.Fs);
            %n(:,[1 end]) = [];
            
            M = mean(mean(n(:,find(dt2 <0)),2));
            SD = std(mean(n(:,find(dt2 <0)),2),0,1);
            M = repmat(M,size(n));
            SD = repmat(SD,size(n));
            n = (n-M)./SD;
            n(isnan(n)) = 0;
            n(isinf(n)) = 0;
            
            if any(isnan(n))
                error('NaNs detected');
            end;
            
            c3=c3+1;
            M3(c3,:) = sum(n,1)./size(n,1);
            
        end;
        
    end;
    
    %%
    x1 = [dat.FRce{:}]';
    x2 = [dat.FRc{:}]';
    x3 = [dat.FRe{:}]';
    
    x4 = dat.nItemsCE;
    x5 = dat.nItemsC;
    x6 = dat.nItemsE;
    
    %%
    if ~isempty(x1)
        FRc = [FRce;x1];
    end;
    
    if ~isempty(x2)
        FRc = [FRc;x2];
    end;
    
    if ~isempty(x3)
        FRe =  [FRe;x3];
    end;
    
    if ~isempty(x4)
        nItemsCE = [nItemsC;x4];
    end;
    
    if ~isempty(x5)
        nItemsC = [nItemsC;x5];
    end;
    
    if ~isempty(x6)
        nItemsE = [nItemsE;x6];
    end;
    
end;
