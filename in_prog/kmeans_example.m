%% load data

load fisheriris;
dat = meas;%(:,1:2);
clear meas;

%% simulate data
% dat = randn(100,2);
% dat(1:25,1) = dat(1:25,1)+3;
% dat(1:25,2) = dat(1:25,2)-4;

%%
nc = 3;
n = size(dat);
niter = 100;
%%
for zt = 1:6
    P = cell(niter,1);
    KM = cell(niter,1);
    
    W = zeros(niter,1);
    f = 0;
    ct = 0;
    while f <1
        
        ct = ct +1;
        
        % initial seed
        if ct ==1
            
            C = zeros(n(1),1);
            for jt = 1:n(1) % randomly assign a label to each observation
                cf = randperm(nc);
                C(jt) = cf(1);
            end;
        end;
        iP= C;%initial partition
        
        % compute the k-means
        c = zeros(nc,n(2));
        for kt = 1:nc% loop over k clusters
            for jt = 1:n(2)% loop over dims
                c(kt,jt) = mean(dat(find(C == kt),jt)); %
            end;
        end;
        KM{ct} = c;
        
        % re-assign cluster labels based on distance from centroid
        C = zeros(n(1),1);
        d3 = zeros(n(1),1);
        for it = 1:n(1)
            
            x = dat(it,:);
            
            d1 = zeros(nc,1);
            for jt = 1:nc
                for nt = 1:length(x)
                    dum(nt) = (x(nt)-c(jt,nt))^2;
                end;
                d1(jt) = sqrt(sum(dum));
            end;
            
            [~,ix] = min(d1);
            
            C(it) = ix;
            
            d3(it) = sum( d1 );
            
        end;
        P{ct} = C;
        
        % compute the variability within cluster
        w = zeros(nc,1);
        for kt = 1:nc
            
            % make pair-wise combinations
            c = 0;
            ix1 = find(C == kt );
            np = length(find(C==kt))*((length(find(C==kt))-1)/2);
            pix = zeros(np,2);
            for it = 1:length(ix1)
                ix2 = setdiff(ix1,ix1(it));
                for jt = 1:length(ix2)
                    if sum(ismember(sort([ix1(it) ix2(jt)],2),sort(pix,2),'rows')) == 0
                        c = c+1;
                        pix(c,:) = [ix1(it) ix2(jt)];
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
            
            % compute the variability within cluster
            w(kt) = 1/length(ix1) * sum(d2);
            
        end;
        
        W(ct) = sum(w);
        
        if (ct-1 >1) && (isequal(P{ct},P{ct-1})) && ( isequal(W(ct),W(ct-1)) )
            break
        end;
    end;
    
    %%
    C = {'y' 'm' 'g'};
    
    figure;
    subplot(2,4,1);
    plot(dat(:,1),dat(:,2),'k.');
    
    subplot(2,4,2);
    hold on;
    for kt = 1:nc
        plot(dat(find(P{1}==kt),1),dat(find(P{1}==kt),2),'.','Color',C{kt});
    end;
    for kt = 1:nc
        plot(KM{1}(kt,1),KM{1}(kt,2),'ko','MarkerFaceColor',C{kt},'MarkerSize',10);
    end;    
    title(W(1));
    
    subplot(2,4,3);
    hold on;
    for kt = 1:nc
        plot(dat(find(P{2}==kt),1),dat(find(P{2}==kt),2),'.','Color',C{kt});
    end;
    for kt = 1:nc
        plot(KM{2}(kt,1),KM{2}(kt,2),'ko','MarkerFaceColor',C{kt},'MarkerSize',10);
    end;  
    title(W(2));
    
    subplot(2,4,4);
    hold on;
    for kt = 1:nc
        plot(dat(find(P{zt}==kt),1),dat(find(P{zt}==kt),2),'.','Color',C{kt});
    end;
    for kt = 1:nc
        plot(KM{ct}(kt,1),KM{ct}(kt,2),'ko','MarkerFaceColor',C{kt},'MarkerSize',10);
    end;      
    title(W(ct));
    
    subplot(2,4,5:6);
    plot([1:ct], W(1:ct) , 'ks-');
    ylabel('Squared euclidean distance');
    xlabel('Iteration');
end;