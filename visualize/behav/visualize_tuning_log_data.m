%%
f.p2logf = '/home/rouxf/Data/EM/Log/P01/';
f.logf = 'P01_24-May-2015_00_00_0000000_log_tune.txt';
p.ncols = 8;

[LogDat] = getNewLogDataTune(f,p)

stim_id = str2double(LogDat.dat(2:end,1));

sel_idxf = find(strcmp(LogDat.dat(2:end,3),'f')==1);
sel_idxp = find(strcmp(LogDat.dat(2:end,3),'p')==1);
sel_idxa = find(strcmp(LogDat.dat(2:end,3),'a')==1);

[n,x] = hist(stim_id,sort(unique([stim_id])));
figure;
subplot(221);
plot(x,n,'b.');axis tight;
ylim([0 max(n)+2]);
xlabel('Stimulus ID');
ylabel('Count');
set(gca,'XTick',[1 max(x)]);
set(gca,'YTick',[0 max(n)]);
subplot(222);
nf = length(sel_idxf);
np = length(sel_idxp);
na = length(sel_idxa);
hold on;
bar(1,nf,'FaceColor',[.9 .75 0]);
bar(2,np,'FaceColor',[.9 .0 0]);
bar(3,na,'FaceColor',[.75 .75 .75]);
ylabel('Count');
set(gca,'XTick',[1 2 3]);
set(gca,'XTickLabel',{'Faces' 'Places' 'Animals'});
subplot(223);
RTs = str2double(LogDat.dat(2:end,end));
hold on;
plot([1],[mean(RTs(sel_idxf))],'s','Color',[.9 .75 0],'MarkerFaceColor',[.9 .75 0]);
plot([2],[mean(RTs(sel_idxp))],'s','Color',[.9 0 0],'MarkerFaceColor',[.9 0 0]);
plot([3],[mean(RTs(sel_idxa))],'s','Color',[.75 .75 .75],'MarkerFaceColor',[.75 .75 .75]);
errorbar(1,mean(RTs(sel_idxf)),std(RTs(sel_idxf)),'Color',[.9 .75 0]);
errorbar(2,mean(RTs(sel_idxp)),std(RTs(sel_idxp)),'Color',[.9 .0 0]);
errorbar(3,mean(RTs(sel_idxa)),std(RTs(sel_idxa)),'Color',[.75 .75 .75]);
ylabel('Reaction times [s]');
set(gca,'XTick',[1 2 3]);
set(gca,'XTickLabel',{'Faces' 'Places' 'Animals'});
subplot(224);
t=str2double(LogDat.dat(2:end,[5:7]));
hold on;
plot([1],[mean(t(:,2)-t(:,1))],'s','Color',[0 0 .9]);
plot([2],[mean(t(:,3)-t(:,2))],'s','Color',[0 0 .9]);
errorbar(1,mean(t(:,2)-t(:,1)),std(t(:,2)-t(:,1)),'Color',[0 0 .9]);
errorbar(2,mean(t(:,3)-t(:,2)),std(t(:,3)-t(:,2)),'Color',[0 0 .9]);
ylabel('Times [s]');
set(gca,'XTick',[1 2]);
set(gca,'XTickLabel',{'\Deltat1' '\Deltat2'});