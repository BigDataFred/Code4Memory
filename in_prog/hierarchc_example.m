%% load data

%load fisheriris meas;
%dat = meas(:,1:2);
%clear meas;

%% simulate data
dat = randn(100,2);
dat(1:25,1) = dat(1:25,1)+3;
dat(1:25,2) = dat(1:25,2)-4;

%%
n = size(dat);

%%
ct = 0;
lin = zeros(n(1)-1,3);
for zt = n(1):-1:2
    
    if zt == n(1)
        cID = 1:n(1);
        cIDo = cID;
    end;
    
    cL = unique(cID);
    
    % make pair-wise inter-cluster pairs
    ix = {};
    for kt = 1:length(cL)
        
        ix{kt} = find(cID == cL(kt));
        
    end;
    
    c = 0;
    np = length(cL)*(length(cL)-1)/2;
    pix = zeros(np,2);
    
    for kt = 1:length(ix)
        
        ix1 = ix{kt};        
        ix2 = setdiff(1:length(ix),kt);
        
        for it = 1:length(ix2)
            ix3 = ix{ix2(it)};
            for jt = 1:length(ix1)
                for nt = 1:length(ix3)
                    if sum(ismember(sort([ix1(jt) ix3(nt)],2),sort(pix,2),'rows')) == 0
                        c = c +1;
                        pix(c,:) = [ix1(jt) ix3(nt)];
                    end;                    
                end;
            end;
        end;
    end;        
    
    % compute euclidean distance for each observation
    d2 = zeros(size(pix,1),1);
    for it = 1:size(pix,1)
        
        x = zeros(size(dat,2),1);
        y = zeros(size(dat,2),1);
        for jt =1:n(2)
            x(jt) = dat(pix(it,1),jt);
            y(jt) = dat(pix(it,2),jt);
        end;
        
        d2(it) = sqrt(sum((x-y).^2));
        
    end;
    
    [v,ix] = min(d2);
    
    ct=ct+1;
    cID(pix(ix,:)) = max(cIDo)+ct;
    
    lin(ct,:) = [pix(ix,:) v];
    
end;
cID - max(cIDo)