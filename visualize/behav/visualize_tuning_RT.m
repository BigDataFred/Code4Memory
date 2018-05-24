function visualize_tuning_RT(bdat)
figure;
subplot(131);
hold on;
bar(bdat.RT(bdat.s_idx));
plot(bdat.RT(bdat.s_idx),'r.');
plot(find(ismember(bdat.s_idx,bdat.oL)),bdat.RT(bdat.s_idx(ismember(bdat.s_idx,bdat.oL))),'go');
axis tight;

subplot(132);
[n,x] = hist(log10(bdat.RT(bdat.s_idx(~ismember(bdat.s_idx,bdat.oL)))),linspace(min(log10(bdat.RT(bdat.s_idx(~ismember(bdat.s_idx,bdat.oL))))),max(log10(bdat.RT(~ismember(bdat.s_idx,bdat.oL)))),100));
bar(x,n);
axis tight;

subplot(133);
sel_idx = setdiff(1:length(bdat.RT),bdat.oL);
idx1 = find(strcmp(bdat.cond,'a'));
idx2 = find(strcmp(bdat.cond,'p'));
idx3 = find(strcmp(bdat.cond,'f'));
y1 = bdat.RT(intersect(idx1,sel_idx));
y2 = bdat.RT(intersect(idx2,sel_idx));
y3 = bdat.RT(intersect(idx3,sel_idx));
hold on;
h(1) = bar(1,[mean(y2)]);
h(2) = bar(2,[mean(y3)]);
plot([1 1],[mean(y2) mean(y2)+std(y2)/sqrt(length(y2)-1)],'k');
plot([2 2],[mean(y3) mean(y3)+std(y3)/sqrt(length(y3)-1)],'k');
plot([.9 1.1],[mean(y2)+std(y2)/sqrt(length(y2)-1) mean(y2)+std(y2)/sqrt(length(y2)-1)],'k');
plot([1.9 2.1],[mean(y3)+std(y3)/sqrt(length(y3)-1) mean(y3)+std(y3)/sqrt(length(y3)-1)],'k');
xlim([0 3]);
ylim([1 max([mean(y2)+std(y2);mean(y3)+std(y3)])]);
set(h(1),'FaceColor','k','EdgeColor','k');
set(h(2),'FaceColor','w','EdgeColor','k');