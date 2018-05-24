function launch_preproc_chan(chan_idx)
%%
addpath('~rouxf/MATLAB/toolboxes/fieldtrip-20160309/');
ft_defaults;
addpath('~rouxf/MATLAB/utils/logfile_readers/');    
%%
p2d = '/home/rouxf/Data/EM/Neuralynx_data/tune/P01/2016-05-25_10-12-44/';
%%
[hdr]   = ft_read_header(p2d);
[event] = ft_read_event(p2d);
%%
f.path2file = '/home/rouxf/Data/EM/Log/P01/';
f.filename = 'P01_25-May-2016_10_13_3612100_log_tune.txt';
p.ncols = 8;
[LogDat] = getNewLogDataTune(f,p);
%[LogDat] = getNewLogDataTune();
load([f.path2file,'P01_tuning_param_25-May-2016_10:13:44.mat']);
%%
make_figure = 'no';
%%
tc = zeros(length(params.idx),1);
for it = 1:length(params.idx)
    
   tc(it) = params.tc(params.idx(params.idx(it)));
      
end;

ntrl = length(LogDat.dat(2:end,1));
c = zeros(length(params.tc),1);
idx2 = zeros(ntrl,1);
for it = 1:ntrl
    
    idx = find([event(:).value] == tc(it));
    
    c(find(tc(it)==params.tc)) = c(find(tc(it)==params.tc)) + 1;
    
    idx2(it) = idx(c(find(tc(it)==params.tc)));
    
end;
ttl_idx = idx2;

dx = diff([tc(1:ntrl) [event(ttl_idx).value]'],[],2);
if sum(dx) ~=0
    error('trigger assignment must match');
end;
%%
cond = LogDat.dat(2:end,3);
stimID = str2double(LogDat.dat(2:end,1));
resp = str2double(LogDat.dat(2:end,4));
RT = str2double(LogDat.dat(2:end,end));
%%
M = nanmean(RT);
SD = nanstd(RT);
z = (RT-M)./SD;
oL = find(z >= 3.5);

[v,s_idx] = sort(RT);
%%
if strcmp(make_figure,'yes')
figure;
subplot(131);
hold on;
bar(v);
plot(v,'r.');
plot(find(ismember(s_idx,oL)),v(ismember(s_idx,oL)),'go');
axis tight;

subplot(132);
[n,x] = hist(log10(v(~ismember(s_idx,oL))),linspace(min(log10(v(~ismember(s_idx,oL)))),max(log10(v(~ismember(s_idx,oL)))),100));
bar(x,n);
axis tight;

subplot(133);
sel_idx = setdiff(1:length(RT),oL);
idx1 = find(strcmp(cond,'a'));
idx2 = find(strcmp(cond,'p'));
idx3 = find(strcmp(cond,'f'));
y1 = RT(intersect(idx1,sel_idx));
y2 = RT(intersect(idx2,sel_idx));
y3 = RT(intersect(idx3,sel_idx));
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
end;
%%
selChan = [ 1:48 ];

fn = cell(length(selChan),1);
for it = 1:length(selChan)
    fn{it} = dir([p2d,'A/A',num2str(selChan(it)),'.Ncs']);
end;

x = double([event(ttl_idx).sample]);
x = x';
%%

trl             = [x-1*hdr.Fs x+3*hdr.Fs -1*hdr.Fs*ones(size(x,1),1)];
%%
CSC_preproc = cell(length(chan_idx),1);
for it = 1:length(chan_idx)
	dataset = [p2d,'A/',fn{chan_idx(it)}.name];
	[CSC_preproc{it}] = preproc_LFP(dataset,trl);
end;
[CSC_preproc] = ft_appenddata([],CSC_preproc{:});
%%
savepath = [p2d,'A/preproc/'];
savename = ['A',num2str(chan_idx(1)),':',num2str(chan_idx(end)),'_preproc.mat'];

save([savepath,savename],'CSC_preproc');

