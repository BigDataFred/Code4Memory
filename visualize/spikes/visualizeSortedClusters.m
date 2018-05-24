%%
[pc,score,latent,tsquare] = princomp(dum.wavf);

cID = unique(dum.assignedCluster);
selIx = {};
for kt = 1:length(cID)
    selIx{kt} = find(dum.assignedCluster == cID(kt));
end;
c = {'r' 'b' 'c' 'm' 'y' 'g'};

figure;
subplot(211);
plot(score(:,1),score(:,2),'ko');
subplot(212);
plot(score(:,1),score(:,2),'ko');
hold on;
for kt = 1:length( selIx )
    plot(score(selIx{kt},1),score(selIx{kt},2),'.','Color',c{kt});
end;

figure;
subplot(121);
hold on;
for kt = 1:length( selIx )
    plot(linspace(0,2,64),dum.wavf(selIx{kt},:),'Color',c{kt});
end;
subplot(122);
hold on;
for kt = 1:length( selIx )
    plot(linspace(0,2,64),mean(dum.wavf(selIx{kt},:),1),'Color',c{kt},'LineWidth',3);
end;

p = [];
cnt = 0;
for it = 1:length(cID)
    for jt = it+1:length(cID)
      cnt = cnt+1;
      p(cnt,:) = [it jt];
    end;
end;

figure;
m = 3;
n = ceil(size(p,1)/m);
for it = 1:size(p,1)
    spikes1 = (dum.wavf_decorr(selIx{p(it,1)},:));
    spikes2 = (dum.wavf_decorr(selIx{p(it,2)},:));
    
    [m1,m2, residuals1, residuals2, overlap,d ] = projectionTest( spikes1,spikes2 );
    
    [bc1,fhat1,h1]=estimatePDF(residuals1,30);
    [bc2,fhat2,h2]=estimatePDF(residuals2,30);
    
    %this is used for calculating Rsquare
    dist1 = normpdf(bc1,0,1);
    dist2 = normpdf(bc2-d,0,1);
    [Rsquare1] = calcRsquare( fhat1, dist1 );
    [Rsquare2] = calcRsquare( fhat2, dist2 );
    
    %more sampling points for plotting purposes
    xDist=-4:.1:4;
    distPlot=normpdf(xDist,0,1);
    
    subplot(n,m,it);    
    bar(bc1, fhat1,  c{p(it,1)});
    hold on
    bar(bc2, fhat2, c{p(it,2)});
    h = [];
    h(1) = plot(xDist,distPlot,'Color',c{p(it,1)},'LineWidth',2.5);
    h(2) = plot(xDist+d,distPlot,'Color',c{p(it,2)},'LineWidth',2.5);
    legend(h,num2str(Rsquare1),num2str(Rsquare2));
    hold off
    ylim([0 0.5]);
    title(num2str(d));
end;

%%
p = [];
cnt = 0;
for it = 1:size(dum2,2)
    for jt = it+1:size(dum2,2)
        cnt=cnt+1;
        p(cnt,:) = [it jt];        
    end;
end;

n = size(p,1)/(size(dum2,2)-1);
m = size(p,1)/n;
n = (size(dum2,2)-1);
figure;
cnt = 0;
s = 0;
for it = 1:size(p,1)
    cnt = cnt+1;
    subplot(n,m,cnt);
    hold on;
    %plot(dum2(:,p(it,1)),dum2(:,p(it,2)),'o','Color',[.75 .75 .75]);    
    for jt = 1:length(selIx)
       plot(dum2(selIx{jt},p(it,1)),dum2(selIx{jt},p(it,2)),'.','Color',c{jt});  
    end;
    axis off;axis tight;
    if mod(cnt,9) == 0
        s = s+1;
        cnt = cnt+s;
    end;
end;

%%
T = length(dataSamples)/Fs;
W = 1/T;
TW = T*W;
k = 2*TW-1;

params                  = [];
params.pad              = 0;
params.Fs               = Fs;
params.fpass            = [0 30];
params.tapers           = [TW k];

[S,f] = mtspectrumc( gradient(dataSamples)', params );