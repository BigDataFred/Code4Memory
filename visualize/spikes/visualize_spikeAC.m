function visualize_spikeAC(AC)

h = bar(AC.time{1},AC.avg,'LineWidth',1);

set(h,'FaceColor',[.75 .75 .75],'EdgeColor','k');
axis tight;
axis(gca,'tight');
ylabel(gca,'Conditional rate (spikes/s)');
xlabel(gca,'\Deltat (s)');

return;