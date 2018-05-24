p2f = '/home/rouxf/res/tuning/P04/2016-10-17_16-37-42/';
files = dir([p2f,'MUA_data_CSC_L*.ncs_stimlocked.mat']);
%%
ch = zeros(1,length(files)); 
for it = 1:length(files) 
    id = [];
    id(1) = max(regexp(files(it).name,'CSC_'));
    id(2) = min(regexp(files(it).name,'.ncs'));
    
    ch(it) =  str2double(  files(it).name(id(1)+5:id(2)-1) );
end;

[ch,s_idx] = sort(ch);

files  = files(s_idx);

%%

dat = cell(1,length(files));

for it = 1:length(files)
    
    [dat{it}] =load([p2f,files(it).name]);
    
end;
%%
chanID = ch;
%[chanID] = makeCSClabels();
n = length(chanID)/8;
%%
t_idx = find(dat{1}.mua.time >= -.35 & dat{1}.mua.time <0);

wsz = 51;
halfwin = floor(wsz/2);
wst = 3;
gw = normpdf(linspace(-halfwin,halfwin,wsz),0,wst);
gw = gw./sum(gw);

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
    y = dat{it}.mua.avg;
    M = mean(y(t_idx));
    SD = std(y(t_idx));    
    z = (y-M)./(SD);
    zM(it) = max(abs(z));
    %y = z;
    y=[];
    y = squeeze(dat{it}.mua.trial);
    
    for jt = 1:size(y,1)
        dum = conv([y(jt,:) y(jt,:) y(jt,:)],gw,'same');
        y(jt,:) = dum(size(y,2)+1:size(y,2)+size(y,2));
    end;
       
    y = y - median(mean(y(:,t_idx),2),1);
    y = mean(y,1);
    
    lm = [lm min(y) max(y)];
    
    plot(dat{it}.mua.time,y);    
    
    axis tight;
    xlim([-.5 1]);
    xlabel('Time (s)');
    ylabel('MUA (a.u.)');
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
set(a,'YLim',[min(lm) max(lm)]);


