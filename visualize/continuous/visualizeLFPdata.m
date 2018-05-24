if isempty(gcp('nocreate'))
    parpool(32);
end;
%%
p2f = '/home/rouxf/res/tuning/P02/2016-07-09_11-13-28/';
files = dir([p2f,'LFP_data_CSC_*.ncs_stimlocked.mat']);
%%
f = [];
f.p2logf = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/P02/log/Tuning/'];
f.logf = 'P02_TS01_log_ctune_09-Jul-2016_11_40_xx.txt';
p.ncols = 8;

[LogDat] = getNewLogDataTune(f,p);
%%
ch = {};%zeros(1,length(files)); 
parfor it = 1:length(files) 
    id = [];
    id(1) = max(regexp(files(it).name,'CSC_'));
    id(2) = min(regexp(files(it).name,'.ncs'));
    
    ch{it} =  str2double(  files(it).name(id(1)+5:id(2)-1) );
    if isnan(ch{it})
        ch{it} =  files(it).name(id(1)+4:id(2)-1);
    end;
end;

[ch,s_idx] = sort(ch);

files  = files(s_idx);

%%
dat = cell(1,length(files));

parfor it = 1:length(files)
    
    [dat{it}] =load([p2f,files(it).name]);
    
end;
%%
chanID = ch;
%[chanID] = makeCSClabels();
n = length(chanID)/8;
%%
t_idx = find(dat{1}.tlck.time >= -.5 & dat{1}.tlck.time <0);

figure(1);
c = 0;
a = zeros(1,length(dat));
lm = [];
zM = zeros(1,length(dat));
for it = 1:length(dat)
    c = c+1;
    subplot(n,8,c);
    a(it) = gca;
    hold on;
    y = dat{it}.tlck.avg;
    M = mean(y(t_idx));
    SD = std(y(t_idx));    
    z = (y-M)./(SD);
    zM(it) = max(abs(z));
    %y = z;
    lm = [lm min(y) max(y)];
    
    plot(dat{it}.tlck.time,y);    
    
    axis tight;
    xlim([-.5 1]);
    xlabel('Time (s)');
    ylabel('Amp. (\muV)');
    title([chanID(it)]);
end;
for it =1:length(a)
    plot(a(it),[0 0],[min(lm) max(lm)],'r');
end;

zM = lm(find(sign(lm)==1));
zM = (zM -mean(zM))./std(zM);

sel_idx = find(zM >1.5);
for it = 1:length(sel_idx)
    area(a(sel_idx(it)),[-.5 1],[max(lm) max(lm)],min(lm),'FaceColor',[.9 0 0],'FaceAlpha',.1);
end;
set(gcf,'Color','w');
%set(a,'YLim',[min(lm) max(lm)]);

% figure;
% hold on;
% plot(sort(zM),'wo','MarkerFaceColor','b')
% plot([1 length(zM)],[1.5 1.5],'r--');
% xlabel('Channel #');
% ylabel('z-score [\sigma]');
% set(gcf,'Color','w');

%%
[~,s_idx] = sort(LogDat.RT);
    
figure(3);
c = 0;
for it = 1:length(sel_idx)
    c = c+1;
    subplot(ceil(length(sel_idx)/2),2,c);
    hold on;
    imagesc(dat{sel_idx(it)}.tlck.time,1:size(dat{sel_idx(it)}.tlck.trial,1),squeeze(dat{sel_idx(it)}.tlck.trial(s_idx,:,:)));
    plot([0 0],[1 size(dat{sel_idx(it)}.tlck.trial,1)],'w');
    axis tight;
    ca = caxis;
    caxis(ca./10);
    axis xy;
    xlim([-.5 1]);
    xlabel('Time (s)');
    ylabel('Trial #');
    title(chanID(sel_idx(it)));
end;
set(gcf,'Color','w');
%%
figure(3);
c = 0;
ca = zeros(length(dat),2);
a = zeros(1,length(dat));
zM = zeros(1,length(dat));
zm = zeros(1,length(dat));
for it = 1:length(dat)
    
    Y = squeeze(mean(dat{it}.powH.powspctrm,1));
    zm(it) = max(max(Y));
    t_idx = find(dat{it}.powH.time >=-.35 & dat{it}.powH.time <0);
    M = repmat(mean(Y(:,t_idx),2),[1 size(Y,2)]);
    Y = (Y - M)./M;           
    zM(it) = round(max(max(Y))*10)/10;
    
    c = c+1;
    subplot(n,8,c);
    a(it) = gca;
    hold on;
    imagesc(dat{it}.powH.time,dat{it}.powH.freq,Y);
    text(.5,100,num2str(zM(it)),'Color','w');
    axis xy;
    ca(it,:) = caxis;
    plot([0 0],[min(dat{it}.powH.freq) max(dat{it}.powH.freq)],'w');
    axis tight;
    xlim([-.5 1]);
    
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title([chanID(it)]);
end;
zM = (zM - mean(zM))./std(zM);
sel_idx = find(zM > 2);

for it = 1:length(sel_idx)
    box(a(sel_idx(it)),'on');
    set(a(sel_idx(it)),'LineWidth',3);
    set(get(a(sel_idx(it)),'Xaxis'),'Color',[.9 .5 0]);
    set(get(a(sel_idx(it)),'Yaxis'),'Color',[.9 .5 0]);
end;

for it = 1:length(a)
    caxis(a(it),[min(min(ca)) max(max(ca))]./10);
end;
set(gcf,'Color','w');
%%
figure(4);
c = 0;
ca = zeros(length(dat),2);
a = zeros(1,length(dat));
for it = 1:length(dat)
    c = c+1;
    
    Y = squeeze(mean(dat{it}.powL.powspctrm,1));
    t_idx = find(dat{it}.powL.time >=-.35 & dat{it}.powL.time <0);
    M = repmat(mean(Y(:,t_idx),2),[1 size(Y,2)]);
    Y = (Y - M)./M;
    zM(it) = round(max(max(Y))*10)/10;
    
    subplot(n,8,c);
    a(it) = gca;
    hold on;
    imagesc(dat{it}.powL.time,dat{it}.powL.freq,Y);
    text(.5,10,num2str(zM(it)),'Color','w');
    ca(it,:) = caxis;
    axis xy;
    plot([0 0],[min(dat{it}.powL.freq) max(dat{it}.powL.freq)],'w');
    axis tight;
    xlim([-.5 1]);
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title([chanID(it)]);
end;

sel_idx = find(zM > 2);

for it = 1:length(sel_idx)
    box(a(sel_idx(it)),'on');
    set(a(sel_idx(it)),'LineWidth',3);
    set(get(a(sel_idx(it)),'Xaxis'),'Color',[0 .9 0]);
    set(get(a(sel_idx(it)),'Yaxis'),'Color',[0 .9 0]);
end;

for it = 1:length(a)
    caxis(a(it),[min(min(ca)) max(max(ca))]./10);
end;
set(gcf,'Color','w');