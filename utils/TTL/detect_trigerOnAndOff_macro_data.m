%%
tAx = trigDat.time{1};
y = trigDat.trial{1};
dy = diff(y);

sel_idx = find(tAx >=5 & tAx <=15);

onIX = find(sign(dy(sel_idx))==-1)+1;
offIX = find(sign(dy(sel_idx))==1);

figure;
subplot(211);
a =gca;
hold on;
plot(tAx(sel_idx),y(sel_idx),'b.-');
plot(tAx(sel_idx(onIX)),y(sel_idx(onIX)),'rs');
plot(tAx(sel_idx(offIX)),y(sel_idx(offIX)),'bs');
axis tight;

subplot(212);
a = [a gca];
plot(tAx(sel_idx(2:end)),dy(sel_idx(1:end-1)),'b.-');
axis tight;

set(a,'Xlim',[8.35 8.9]);