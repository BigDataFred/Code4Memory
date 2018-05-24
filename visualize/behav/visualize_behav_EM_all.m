%%
addpath(genpath('/media/rouxf/rds-share/Fred/code/mcode/custom/'));

%%
acc = {};
nBlocks= {};
RTs = {};

pID = {'P2B';'P4B';'P5B';'P7B';'P22A';'P23A'};

%%
rpath = '/media/rouxf/rds-share/iEEG_DATA/MICRO/P02/fvSpEM/';
[RTs{1},acc{1},nBlocks{1}] = extract_behavior_EM(rpath);

%%
rpath = '/media/rouxf/rds-share/iEEG_DATA/MICRO/P04/fvSpEM/';
[RTs{2},acc{2},nBlocks{2}] = extract_behavior_EM(rpath);

%%
rpath = '/media/rouxf/rds-share/iEEG_DATA/MICRO/P05/fvSpEM/';
[RTs{3},acc{3},nBlocks{3}] = extract_behavior_EM(rpath);

%%
rpath = '/media/rouxf/rds-share/iEEG_DATA/MICRO/P22AMS/cnEM/';
[RTs{4},acc{4},nBlocks{4}] = extract_behavior_EM(rpath);

%%
rpath = '/media/rouxf/rds-share/iEEG_DATA/MICRO/P23AMS/cnEM/';
[RTs{5},acc{5},nBlocks{5}] = extract_behavior_EM(rpath);

%%
rpath = '/media/rouxf/rds-share/iEEG_DATA/MICRO/P23AMS/cnEM/';
[RTs{5},acc{5},nBlocks{5}] = extract_behavior_EM(rpath);

%%
Macc= [];
SEMacc = [];
d = [];
nB = [];
rt = [];
SEMrt = [];

for it = 1:length(acc)
    
    Macc(it,:) = mean(acc{it},1);
    SEMacc(it,:) = std(acc{it},0,1)/sqrt(size(acc{it},1)-1);
    
    d(it,:)  = [median([nBlocks{it}{2,:}]) std([nBlocks{it}{2,:}])];
    nB(it,:) = [median([nBlocks{it}{1,:}]) std([nBlocks{it}{1,:}])];
    
    rt(it,:) = [mean(RTs{it}{1}) mean(RTs{it}{2}) mean(RTs{it}{3})];
    SEMrt(it,:) = [std(RTs{it}{1})/sqrt(length(RTs{it}{1})-1) std(RTs{it}{2})/sqrt(length(RTs{it}{2})-1) std(RTs{it}{3})/sqrt(length(RTs{it}{3})-1)];
    
end;
SEMrt(isnan(SEMrt)) = 0;

[~,six] = sort(Macc(:,1));

figure;
subplot(221);
a = gca;
hold on;
errorbar(1:size(Macc,1),Macc(six,1),SEMacc(six,1),'bs-','LineWidth',1);
errorbar(1:size(Macc,1),mean(Macc(six,2:3),2),mean(SEMacc(six,2:3),2),'rs-','LineWidth',1);
plot([1 size(Macc,1)],[.5 .5],'k--');
subplot(222);
a = [a gca];
hold on;
errorbar(rt(six,1),SEMrt(six,1),'bs-','LineWidth',1);
errorbar(nanmean(rt(six,2:3),2),nanmean(SEMrt(six,2:3),2),'rs-','LineWidth',1);
subplot(223);
a = [a gca];
errorbar(nB(six,1),nB(six,2),'bs-','LineWidth',1);
subplot(224);
a = [a gca];
errorbar(d(six,1),d(six,2),'bs-','LineWidth',1);

axis(a,'tight');
set(a(1),'YLim',[0 1]);
set(a(1),'YTick',[.25 .5 .75]);
set(a,'XLim',[0 size(Macc,1)+1]);
set(a,'XTick',1:size(Macc,1));
set(a,'XTickLabel',pID(six));

title(a(1),'Accuracy');
title(a(2),'Reaction times');
title(a(3),'Blocks/session');
title(a(4),'Pairs/block');

ylabel(a(1),'[%]');
ylabel(a(2),'[s]');
ylabel(a(3),'Count');
ylabel(a(4),'Count');
