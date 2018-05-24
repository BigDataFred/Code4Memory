%%
SE1 = nanstd(M1,0,1)./sqrt(size(M1,1)-1);
Y1 = nanmean(M1,1);
SE2 = nanstd(M2,0,1)./sqrt(size(M2,1)-1);
Y2 = nanmean(M2,1);
SE3 = nanstd(M3,0,1)./sqrt(size(M3,1)-1);
Y3 = nanmean(M3,1);

h = [];
figure;
a =gca;
hold on;
plot([0 0],[min(Y1-SE1) max(Y1+SE1)],'k');
plot([2000 2000],[min(Y1-SE1) max(Y1+SE1)],'k');
%errorbar(dt2(1:end-1)+250,Y1,SE1,'ms-','MarkerFaceColor','m');
%errorbar(dt2(1:end-1)+250,Y2,SE2,'gs-','MarkerFaceColor','g');
jbfill(dt2(1:end-1)+250,Y1+SE1,Y1-SE1,[255 127 80]./255,[255 127 80]./255,0,.25);
hold on;
h(1) = plot(dt2(1:end-1)+250,Y1,'Color',[255 127 80]./255);

jbfill(dt2(1:end-1)+250,Y2+SE2,Y2-SE2,[46 139 87]./255,[46 139 87]./255,0,.25);
h(2) = plot(dt2(1:end-1)+250,Y2,'Color',[46 139 87]./255);

jbfill(dt2(1:end-1)+250,Y3+SE3,Y3-SE3,[255 0 0]./255,[255 127 80]./255,0,.25);
h(3) = plot(dt2(1:end-1)+250,Y3,'Color',[255 0 0]./255);
axis tight;

for it = 1:length(a)
    xlabel(a(it),'Time [ms]');
    ylabel(a(it),'Z-score [\sigma]');
end;

set(gcf,'Color','w');


%%

figure;
subplot(131);
a =gca;
hold on;
for it = 1:size(nItemsCE,2)
    M  = mean( nItemsCE(:,it));
    SE = std(nItemsCE(:,it))/sqrt(length(nItemsCE(:,it))-1);
        
    plot([it-.2 it+.2],[M M],'r','LineWidth',3);
    plot([it-.2 it+.2],[M-SE M-SE],'b','LineWidth',3);
    plot([it-.2 it+.2],[M+SE M+SE],'b','LineWidth',3);
    plot([it it]-.2,[M-SE M+SE],'b','LineWidth',3);
    plot([it it]+.2,[M-SE M+SE],'b','LineWidth',3);
    plot( it*ones(1,size(nItemsCE(:,it),1)) , nItemsCE(:,it),'o','Color',[46 139 87]./255,'MarkerFaceColor','k');
end;
for it = 1:size(nItemsCE,1)
    plot([1 2],nItemsCE(it,:),'k-')
end;
xlim([0 size(nItemsCE,2)+1]);
title('CE-Units');

subplot(132);
a =[a gca];
hold on;
for it = 1:size(nItemsC,2)
    M  = mean( nItemsC(:,it));
    SE = std(nItemsC(:,it))/sqrt(length(nItemsC(:,it))-1);
        
    plot([it-.2 it+.2],[M M],'r','LineWidth',3);
    plot([it-.2 it+.2],[M-SE M-SE],'b','LineWidth',3);
    plot([it-.2 it+.2],[M+SE M+SE],'b','LineWidth',3);
    plot([it it]-.2,[M-SE M+SE],'b','LineWidth',3);
    plot([it it]+.2,[M-SE M+SE],'b','LineWidth',3);
    plot( it*ones(1,size(nItemsC(:,it),1)) , nItemsC(:,it),'o','Color',[46 139 87]./255,'MarkerFaceColor','k');
end;
for it = 1:size(nItemsC,1)
    plot([1 2],nItemsC(it,:),'k-')
end;
xlim([0 size(nItemsC,2)+1]);
title('C-Units');

subplot(133);
a=[a gca];
hold on;
for it = 1:size(nItemsE,2)
    M  = mean( nItemsE(:,it));
    SE = std(nItemsE(:,it))/sqrt(length(nItemsE(:,it))-1);
        
    plot([it-.2 it+.2],[M M],'r','LineWidth',3);
    plot([it-.2 it+.2],[M-SE M-SE],'b','LineWidth',3);
    plot([it-.2 it+.2],[M+SE M+SE],'b','LineWidth',3);
    plot([it it]-.2,[M-SE M+SE],'b','LineWidth',3);
    plot([it it]+.2,[M-SE M+SE],'b','LineWidth',3);
    plot( it*ones(1,size(nItemsE(:,it),1)) , nItemsE(:,it),'o','Color',[255 127 80]./255,'MarkerFaceColor','k');
end;
for it = 1:size(nItemsE,1)
    plot([1 2],nItemsE(it,:),'k-')
end;
xlim([0 size(nItemsE,2)+1]);
title('E-Units');

set(a,'XTick',[1 2]);
set(a,'XTickLabel',{'Low' 'High'});

for it = 1:length(a)
    ylabel(a(it),'Average item recall');
end;
set(gcf,'Color','w');


