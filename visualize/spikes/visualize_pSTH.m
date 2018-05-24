function visualize_pSTH(dt,y,C)

%plot(dt,y/(length(y)*(dt(2)-dt(1))),'Color',C);
plot(dt,y,'Color',C);

axis(gca,'tight');
ylabel(gca,'Firing rate (spikes/s)');
xlabel(gca,'Time (s)');