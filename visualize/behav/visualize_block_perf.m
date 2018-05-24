function visualize_block_perf
%%
addpath('/media/rouxf/rds-share/Fred/code/mcode/custom/helper/logfile_readers/');
addpath('/media/rouxf/rds-share/Fred/code/mcode/custom/visualize/behav/');

params                  = [];
params.p                = '/media/rouxf/rds-share/Archive/MICRO/P09/log/EM/';
params.fn               = 'P09_fVSp_02_2017-Aug-30_9_54_1514400_LogFile_EMtask.txt';
params.ntrl             = 99;
params.ncol             = 12;

plot_block_perf(params)

function plot_block_perf(params)
%%
[LogDat] = getNewLogDataEM(params,'RET');

%%
if ~(length(unique(LogDat.log(:,1))) == length(LogDat.log(:,1)))
    error('non unique cue detected');
end;
if ~(length(unique(LogDat.log(:,2))) == length(LogDat.log(:,2)))
    error('non unique stim1 detected');
end;
if ~(length(unique(LogDat.log(:,3))) == length(LogDat.log(:,3)))
    error('non unique stim2 detected');
end;

%%
x = []; for it = 1:size(LogDat.log,1);x(it) = isempty([LogDat.log{it,:}]);end;

%%
id = params.fn(1:regexp(params.fn,'\d{4}')+3);

%%
nB = size(LogDat.idx,1);
pct = zeros(nB,2);
n = zeros(nB,1);
rt = zeros(nB,6);
rt2 = cell(nB,3);
dt1 = zeros(nB,4);
dt2 = cell(nB,2);

%%
for it = 1:nB
    
    if size(LogDat.log,2) == 11
        x1 = str2double( LogDat.log( LogDat.idx(it,1):LogDat.idx(it,2),5:6 ) );
        x2 = (str2double( LogDat.log( LogDat.idx(it,1):LogDat.idx(it,2),7:8) ).*1000);
        x3 = str2double( LogDat.log( LogDat.idx(it,1):LogDat.idx(it,2),9:11) ).*1000;
        
        rt(it,1) = nanmean(x2(:,1));
        rt(it,2) = nanstd(x2(:,1));
        rt(it,3) = nanmean(x2(:,2));
        rt(it,4) = nanstd(x2(:,2));        
        
        rt2{it,1} = x2(:,1);
        rt2{it,2} = x2(:,2);
        
    else
        x1 = str2double( LogDat.log( LogDat.idx(it,1):LogDat.idx(it,2),5:6 ) );
        x2 = (str2double( LogDat.log( LogDat.idx(it,1):LogDat.idx(it,2),7:9) ).*1000);
        x3 = str2double( LogDat.log( LogDat.idx(it,1):LogDat.idx(it,2),10:12) ).*1000;
        
        [i,j] = find(x2 == 1e9);x2(i,j) = NaN;
        
        rt(it,1) = nanmean(x2(:,1));
        rt(it,2) = 3*(nanstd(x2(:,1))/(size(x2,1)-1));
        
        rt(it,3) = nanmean(x2(:,2));
        rt(it,4) = 3*(nanstd(x2(:,2))/(size(x2,1)-1));
        
        rt(it,5) = nanmean(x2(:,3));
        rt(it,6) = 3*(nanstd(x2(:,3))/(size(x2,1)-1));
        
        rt2{it,1} = x2(:,1);
        rt2{it,2} = x2(:,2);
        rt2{it,3} = x2(:,3);
        
    end;
    n(it) = length( LogDat.idx(it,1):LogDat.idx(it,2));
    
    pct(it,1) = nansum(x1(:,1))/length(x1(:,1));
    pct(it,2) = nansum(x1(:,2))/length(x1(:,2));    
    
    dt1(it,1) = median(x3(:,2)-x3(:,1));
    dt1(it,2) = std(x3(:,2)-x3(:,1));
    
    dt1(it,3) = median(x3(:,3)-x3(:,2));
    dt1(it,4) = std(x3(:,3)-x3(:,2));
    
    dt2{it,1} = x3(:,2)-x3(:,1);
    dt2{it,2} = x3(:,3)-x3(:,2);
end;

%%
if sum(n) ~= length(find(x~=1))
    error('trial number must match');
end;

%%
figure;
subplot(321);
hold on;
h(1) = area([1 nB],[1 1],.75);
h(2) = area([1 nB],[.75 .75],.65);
h(3) = area([1 nB],[.65 .65],.55);
h(4) = area([1 nB],[.55 .55],.35);
h(5) = area([1 nB],[.35 .35],.15);

set(h(1),'FaceColor',[.9 0 0]);
set(h(2),'FaceColor',[.9 .5 0]);
set(h(3),'FaceColor',[.55 .55 .55]);
set(h(4),'FaceColor',[0 .5 .9]);
set(h(5),'FaceColor',[0 0 .9]);

[A,H1,H2] = plotyy(1:nB,(mean(pct,2)-1/4)/(1-(1/4)),1:nB,n);
set([H1 H2],'LineWidth',3);
axis(A,'tight');

set(H1,'Color','k');
set(H1,'Marker','s');
set(H2,'Color',[.75 .75 .75]);
set(H2,'Marker','s');

set(A(1),'YTick',[.15:.2:1]);

set((A),'YColor','k');

xlabel(A(1),'Number of encoding blocks','Color','k');
ylabel(A(1),'Percent correct [%]','Color','k');
ylabel(A(2),'Number of scenes / block','Color','k');

title({id});

subplot(323);
hold on;
h = [];
h(1) = plot(1:nB,(pct(:,1)-1/4)/(1-(1/4)),'bs-');
h(2) = plot(1:nB,(pct(:,2)-1/4)/(1-(1/4)),'rs-');
ylabel(gca,'Percent correct [%]');
xlabel(gca,'Number of encoding blocks');
ylim([0 1]);xlim([0 nB+1]);
legend(h,'Item #1','Item #2');

subplot(324);
hold on;
%errorbar(1:nB,rt(:,1),rt(:,2),'ks-');
errorbar(1:nB,rt(:,3),rt(:,4),'bs-');
errorbar(1:nB,rt(:,5),rt(:,6),'rs--');

ylabel(gca,'Reaction times [ms]');
xlabel(gca,'Number of encoding blocks');
xlim([0 nB+1]);

subplot(325);
errorbar(1:nB,dt1(:,1),dt1(:,2),'bs-');
ylabel(gca,'Duration of cue interval [ms]');
xlabel(gca,'Number of encoding blocks');
xlim([0 nB+1]);

subplot(326);
errorbar(1:nB,dt1(:,3),dt1(:,4),'bs-');
ylabel(gca,'Duration of response interval');
xlabel(gca,'Number of encoding blocks');
xlim([0 nB+1]);

set(gcf,'Color','w');