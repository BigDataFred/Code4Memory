%%
addpath('~rouxf/prj/Bham/code/mcode/visualize/spikes/');

p2d = '/media/rouxf/rds-share/iEEG_DATA/MICRO/code4MEM/';
fn = dir([p2d,'*_spikeANDwvfData.mat']);
saveParams.savepath = '/media/rouxf/rds-share/iEEG_DATA/MICRO/code4MEM/clusterFigures/';

%%
t = [-2 5]; % time window of interest
t = t.*1e3; % convert s to ms
dt = t(1):250:t(2);% bins used for the isi range

%%
h = figure;
set(h,'visible','off');
cnt = 0;
for it = 1:length(fn)
    fprintf([num2str(it),'/',num2str(length(fn))]);
    dat = load([p2d,fn(it).name]);
    [~,fn2,~] = fileparts(fn(it).name);
    for jt = 1:size(dat.spkRaster,1)
        cnt = cnt+1;
        saveParams.savename = [fn2,'_unitFig',num2str(cnt),'.fig'];
        visualizeUnits4manualSelection(dat.wvfStats(jt),dat.instFR(jt),dat.spkRaster(jt,:),dt,saveParams);               
    end;
    fprintf('\n')
end;

%%
% cnt = 0;
% sel = [];
% cnt = cnt+1;
% [s] = input('Select unit? [y/n]','s');
% if strcmp(s,'y')
%     sel(cnt) = 1;
% else
%     sel(cnt) = 0;
% end;