%% compute coherence spectra lower freqz
T = floor(length(tIx)/Fs);
W = 3;
TW = T*W;
k = 2*TW-1;

paramsC                     = [];
paramsC.Fs                  = Fs;
paramsC.tapers              = [TW k];
paramsC.pad                 = 2;
paramsC.fpass               = [0 30];
paramsC.err                 = [2 0.001];
paramsC.trialave            = 1;

npairs = size(ICs1{4},1)*(size(ICs1{4},1)-1)/2;
p = zeros(npairs,2);
cnt = 0;
for it = 1:size(ICs1{4},1)
    for jt = it+1:size(ICs1{4},1)
        cnt = cnt+1;
        p(cnt,:) = [it jt];
    end;
end;

[C1] = cell( 1, npairs );
[Cerr1] = cell( 1, npairs );
for it = 1:size(p,1)
    fprintf([num2str(it),'/',num2str(size(p,1))]);

    X = ICs1{4}(p(it,1),:);
    Y = ICs1{4}(p(it,2),:);
    X = reshape(X,[length(tIx) length(X)/length(tIx)]);
    Y = reshape(Y,[length(tIx) length(Y)/length(tIx)]);
    
    [C1{it},~,~,~,~,f1,~,~,Cerr1{it}] = coherencyc( X, Y, paramsC );
    fprintf('\n');
end;

%%
pid = unique(p(:,1));
cnt2 = 0;
h = figure;
for it = 1:length(pid)
    ix = find(p(:,1)==pid(it));
    cnt = ((it-1)*length(pid))+(length(pid)-length(ix));
    for jt = 1:length(ix)        
        cnt = cnt+1;
        cnt2 = cnt2+1;
        subplot(length(pid),length(pid),cnt)
        jbfill(f1,Cerr1{cnt2}(1,:),Cerr1{cnt2}(2,:),[0 0 .9],[0 0 .9],1,.5);
        hold on;        
        plot(f1,C1{cnt2},'y');  
        axis tight;
        title(p(cnt2,2));
        ylabel(p(cnt2,1));
    end;
end;

%% compute coherence spectra higher freqz
T = floor(length(tIx)/Fs);
W = 8;
TW = T*W;
k = 2*TW-1;

paramsC                     = [];
paramsC.Fs                  = Fs;
paramsC.tapers              = [TW k];
paramsC.pad                 = 2;
paramsC.fpass               = [20 170];
paramsC.err                 = [2 0.001];
paramsC.trialave            = 1;

npairs = size(ICs2{4},1)*(size(ICs2{4},1)-1)/2;
p = zeros(npairs,2);
cnt = 0;
for it = 1:size(ICs2{4},1)
    for jt = it+1:size(ICs2{4},1)
        cnt = cnt+1;
        p(cnt,:) = [it jt];
    end;
end;

[C2] = cell( 1, npairs );
[Cerr2] = cell( 1, npairs );
for it = 1:size(p,1)
    fprintf([num2str(it),'/',num2str(size(p,1))]);

    X = ICs2{4}(p(it,1),:);
    Y = ICs2{4}(p(it,2),:);
    X = reshape(X,[length(tIx) length(X)/length(tIx)]);
    Y = reshape(Y,[length(tIx) length(Y)/length(tIx)]);
    
    [C2{it},~,~,~,~,f2,~,~,Cerr2{it}] = coherencyc( X, Y, paramsC );
    fprintf('\n');
end;

%%
pid = unique(p(:,1));
cnt2 = 0;
h = figure;
for it = 1:length(pid)
    ix = find(p(:,1)==pid(it));
    cnt = ((it-1)*length(pid))+(length(pid)-length(ix));
    for jt = 1:length(ix)        
        cnt = cnt+1;
        cnt2 = cnt2+1;
        subplot(length(pid),length(pid),cnt)
        jbfill(f2,Cerr2{cnt2}(1,:),Cerr2{cnt2}(2,:),[0 0 .9],[0 0 .9],1,.5);
        hold on;        
        plot(f2,C2{cnt2},'y');  
        axis tight;
        title(p(cnt2,2));
        ylabel(p(cnt2,1));
    end;
end;