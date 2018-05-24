function plotClusterStats(clusterStat)


figure;

n = length(clusterStat);
c = 0;
for it = 1:n
    c = c+1;
    subplot(n,3,c);
    hold on;
    bar(clusterStat{it}.edges1,clusterStat{it}.n1);
    plot(clusterStat{it}.edges1,clusterStat{it}.yGam1,'r');
    axis tight;
    
    c = c+1;
    subplot(n,3,c);
    hold on;
    sel_ix =(length(clusterStat{it}.tvect)-1)/2+2:length(clusterStat{it}.tvect);
    plot(clusterStat{it}.tvect(sel_ix),clusterStat{it}.Cxx(sel_ix));
    axis tight;
    
    c = c+1;
    subplot(n,3,c);
    hold on;
    plot(clusterStat{it}.f,clusterStat{it}.Pxxn);
    axis tight;   
    
end;
return;