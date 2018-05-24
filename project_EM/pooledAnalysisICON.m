%%
p2d = '/media/rouxf/rds-share/iEEG_DATA/MICRO/code4MEM/logDat/';
logfiles = dir([p2d,'*.mat']);
pctH = [];
pctM = [];
for it = 1:length(logfiles)
    
    dat = load([p2d,logfiles(it).name]);
    ntrl = length(dat.LogDat2.log);
    pctH(it) = length(find(sum(str2double(dat.LogDat2.log(:,5:6)),2)==2))/ntrl;
    pctM(it) = length(find(sum(str2double(dat.LogDat2.log(:,5:6)),2)~=2))/ntrl;
    
end;

%%
p2d = '/media/rouxf/rds-share/iEEG_DATA/MICRO/code4MEM/';
mfiles = dir([p2d,'*.mat']);

PId = {};
ExpM = {};
for it = 1:length(mfiles)
    PId(it) = {mfiles(it).name(regexp(mfiles(it).name,'P\d{1,2}\w{0,3}_'):regexp(mfiles(it).name,'_')-1)};
    ExpM(it) ={ mfiles(it).name(regexp(mfiles(it).name,'_\w{2,5}EM_')+1:regexp(mfiles(it).name,'\d{4}')-2)};
end;
PId = PId';
ExpM = ExpM';

selIdx = 1:length(mfiles);%find(strcmp(ExpM,'fvSpEM'));
mfiles = mfiles(selIdx);

%%
cnt= 0; 
instFR = [];
FR = [];
pval1 = [];
pval2 = [];
pctISI = [];
SNR = [];

%%
for it = 1:length( mfiles )        
    
    [dat] = load([p2d,mfiles(it).name],'uSelIx', 'instFR', 'FR','wvfStats','missIdx','hitIdx','trlENC');
    
    fprintf([num2str(it),'/',num2str(length(mfiles))]);
    
    pct = [];
    snr = [];
    cnt2= 0;
    for jt = 1:length(dat.wvfStats)
        cnt2 = cnt2+1;
        pct(cnt2) = dat.wvfStats{jt}.pct;
        snr(cnt2) = dat.wvfStats{jt}.snr;        
    end;
    selIdx = find( pct <3 );
    SNR = [SNR snr(selIdx)];
    pctISI = [pctISI pct(selIdx)];
    
    pval1 = [pval1;dat.uSelIx(selIdx,3)];    
    pval2 = [pval2;dat.uSelIx(selIdx,2)];
    
    t = -2:0.25:5;
    for jt = 1:length(selIdx)   
        
        cnt = cnt+1;
        % transform instantaneous firing rate
        x = dat.instFR{selIdx(jt)}(2:end-1);        
        z = (x-mean(x(t>=-1 & t <0)))./std(x(t>=-1 & t <0));
        instFR(cnt,:) = z;  
        FR(cnt,:) = mean(dat.FR{selIdx(jt)},1);
    end;
        
    fprintf('\n');
  
    
end;

%%
selIdx = find(strcmp(PId,'P23AMS'));
dat = load([p2d,mfiles(selIdx).name],'spkRaster','FR','wvfStats');

%%
M1 = mean(pctH);
M2 = mean(pctM);
SE1 = std(pctH)/sqrt(length(pctH)-1);
SE2 = std(pctM)/sqrt(length(pctM)-1);


figure;
hold on;
bar([1 2],[M1 M2]);
plot([1 1],[M1 M1+SE1],'k');
plot([.9 1.1],[M1+SE1 M1+SE1],'k');
plot([2 2],[M2 M2+SE2],'k');
plot([1.9 2.1],[M2+SE2 M2+SE2],'k');
axis tight;
xlim([0 3]);
ylim([0 1]);
set(gca,'XTick',[1 2]);
set(gca,'XTickLabel',{'Remembered','Forgotten'});

%%

[h,p] = ttest2(dat.FR{1}(find(ismember(dat.missIdx,dat.trlENC)),3),dat.FR{1}(find(ismember(dat.hitIdx,dat.trlENC)),3));

%%
[h1, crit_p, adj_ci_cvrg, adj_p1]=fdr_bh(pval1,0.025);
[h2, crit_p, adj_ci_cvrg, adj_p1]=fdr_bh(pval2,0.025);

%%
p= [];
c = 0;
for it =1:size(instFR,2)
    c = c+1;
    [~, p(c)] = ttest(instFR(:,it),0);
end;

h = p'<(0.05/size(instFR,2));

%%
selIdx =find(h1==1);

M1 = mean(instFR,1);
SE1 = std(instFR,0,1)/sqrt(size(instFR,1)-1);

figure;
hold on;
%h(1) = area([0.5 1.5],ones(1,2)*max(max([M1+SE1;M2+SE2])),min(min([M1-SE1;M2-SE2])))
%h(2) = area([2.5 3.5],ones(1,2)*max(max([M1+SE1;M2+SE2])),min(min([M1-SE1;M2-SE2])))
%plot([0.5 1.5],ones(1,2)*max(max([M1+SE1])),'k');
%plot([2.5 3.5],ones(1,2)*max(max([M1+SE1])),'k');
%set(h,'FaceColor',[.75 .75 .75],'EdgeColor',[.75 .75 .75],'FaceAlpha',.4,'EdgeAlpha',.4);
jbfill(t(2:end-1),M1-SE1,M1+SE1,[.9 0 0],[.9 0 0],0,.25);
hold on;
plot(t(2:end-1),M1,'Color',[.9 0 0]);
plot([t(2) t(end)],[0 0],'k--');

axis tight;
xlim([-.5 4.75]);
xlabel('Time rel. to cue onset (s)');
ylabel('Average norm. firing rate [\sigma]');

%%
selIdx1 =find(h1==1 & h2 ~=1);
selIdx2 =find(h2==1 & h1 ~=1);
selIdx3 = setdiff(1:length(h1),[find(h1==1);find(h2==1)]);

M1 = mean(instFR(selIdx1,:),1);
SE1 = std(instFR(selIdx1,:),0,1)/sqrt(size(instFR,1)-1);
M2 = mean(instFR(selIdx2,:),1);
SE2 = std(instFR(selIdx2,:),0,1)/sqrt(size(instFR,1)-1);
M3 = mean(instFR(selIdx3,:),1);
SE3 = std(instFR(selIdx3,:),0,1)/sqrt(size(instFR,1)-1);

figure;
hold on;
h = [];
h(1) = area([0.5 1.5],ones(1,2)*max(max([M1+SE1;M2+SE2])),min(min([M1-SE1;M2-SE2])));
h(2) = area([2.5 3.5],ones(1,2)*max(max([M1+SE1;M2+SE2])),min(min([M1-SE1;M2-SE2])));
set(h,'FaceColor',[.75 .75 .75],'EdgeColor',[.75 .75 .75],'FaceAlpha',.3,'EdgeAlpha',.3);
jbfill(t(2:end-1),M1-SE1,M1+SE1,[75 183 33]./255,[75 183 33]./255,0,.25);
jbfill(t(2:end-1),M2-SE2,M2+SE2,[117 49 178]./255,[117 49 178]./255,0,.25);
jbfill(t(2:end-1),M3-SE3,M3+SE3,[50 50 50]./255,[50 50 50]./255,0,.25);

hold on;
plot(t(2:end-1),M1,'Color',[75 183 33]./255);
plot(t(2:end-1),M2,'Color',[117 49 178]./255);
%plot(t(2:end-1),M3,'Color',[50 50 50]./255);
plot([t(2) t(end)],[0 0],'k--');

axis tight;
xlim([-.5 4.75]);
xlabel('Time rel. to cue onset (s)');
ylabel('Average norm. firing rate [\sigma]');



%%
[n,x] = hist(pctISI,[0:.1:max(pctISI)]);
[n2,x2] = hist(SNR,[0:.25:max(SNR)]);

figure;
subplot(121);
bar(x,n);
axis tight;
subplot(122);
bar(x2,n2);
axis tight;













