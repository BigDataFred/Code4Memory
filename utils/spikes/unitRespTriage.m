function [spkResp,bfLab,trlENC,trlRET] = unitRespTriage(pID,expMode,sesh)

%%
addpath('~/prj/Bham/code/mcode/utils/spikes/');

%%
if nargin == 0
    pID = 'P07';
    expMode = 'fVSpEM';
    sesh = '2017-04-30_17-35-52';
end;

%%
rootPath = '/media/rouxf/rds-share/iEEG_DATA/MICRO/';
p2d = [rootPath,pID,'/',expMode,'/',sesh,'/'];
spkFn = [pID,'_',expMode,'_',sesh,'_spkDataStimLockedSegmented.mat'];

[~,fName,ext] = fileparts(spkFn);
saveName = [fName,'manuallyMergedUnits',ext];
chck = dir([p2d,saveName]);

%%
if isempty(chck)
    
    [spkDat] = load([p2d,spkFn]);
    trlENC = find(ismember(spkDat.trlSel,spkDat.trlENC));
    trlRET = find(ismember(spkDat.trlSel,spkDat.trlRET));
    
    %%
    [selInfo,selIx] = VisualizeAndScreeningOfClusterData(spkDat);
    
    %%
    [spkResp,bfLab] = mergeSpikeTimesOfDifferentClusters(spkDat.sortedSpikes,spkDat.trlSel,spkDat.chanLab,selInfo,selIx);
    
    %%
    save([p2d,saveName],'spkResp','bfLab','trlENC','trlRET');
    
else
    load([p2d,saveName]);
end;

%% 
% addpath(genpath('~/tbx/chronux_2_11/spectral_analysis/'));
% 
% params                  = [];
% params.Fs               = 1e3; 
% params.pad              = 0;
% params.fpass            = [0 30];
% params.tapers           = [1 1];
% params.trialave         =1;
% 
% dt = -501:5.001e3;
% for it = 1:length(selIx)       
%     
%     spkResp{selIx(it)}.spikeTimes = spkResp{selIx(it)}.spikeTimes(ismember(spkDat.trlSel,spkDat.trlENC)==1);
%     n = zeros(length(spkResp{selIx(it)}.spikeTimes),length(dt)-2);
%     ts = [];
%     
%     for jt = 1:length( spkResp{selIx(it)}.spikeTimes )
%         ts = spkResp{selIx(it)}.spikeTimes{jt};
%         x= hist(ts,dt);
%         x([1 end]) = [];
%         n(jt,:) = x;
%         
%     end;
%     
%     [S,f,R] = mtspectrumpb(n',params);
%     
%     %if R>=1
%         figure;
%         subplot(321);
%         hold on;
%         plot(linspace(0,2,64),spkResp{it}.wvf,'r');
%         plot(linspace(0,2,64),mean(spkResp{it}.wvf,1),'k');        
%         title(R);
%         subplot(322);
%         hold on;
%         plot([0 0],[0 size(n,1)+1],'r');
%         plot([2e3 2e3],[0 size(n,1)+1],'r');
%         dt2 = dt(2:end-1);
%         for jt = 1:size(n,1)
%             x = dt2(n(jt,:)==1);
%             y = jt*ones(1,length(x));
%             x = [x;x];
%             y = [y-.5;y+.5];
%             line(x,y,'Color','k');
%         end;
%         xlim([-500 5e3]);
%         ylim([0 size(n,1)+1]);
%         subplot(323);
%         plot(f,S./R);
%         subplot(324);
%         hold on;
%         r = [];
%         for jt = 1:size(n,1)
%             r(jt,:) = conv(hist(dt2(n(jt,:)==1),dt2(1):250:dt2(end)),gausswin(5),'same')./sum(gausswin(5));
%         end;
%         fr = sum(r,1)./size(r,1)./0.25;
%         bar(dt2(2):250:dt2(end-1),fr(2:end),'k');
%         plot([0 0],[min(fr) max(fr)],'r');
%         plot([2e3 2e3],[min(fr) max(fr)],'r');
%         xlim([-500 5e3]);
%         axis tight;
%         subplot(3,2,5);
%         xc = zeros(1,1001);
%         for jt = 1:size(n,1)
%             xc = xc + xcorr(n(jt,:),500);
%         end;
%         xc(501) = NaN;
%         xc = xc./jt;
%         plot(-500:500,xc);
%         xlim([-500 500]);
%         subplot(3,2,6);
%         dx = diff([spkResp{it}.spikeTimes{:}]);
%         dx(sign(dx)==-1) = [];
%         [isi,x] = hist(dx,[0:251]);
%         bar(x(1:end-1),isi(1:end-1));
%         xlim([-10 250]);
%         
%         set(gcf,'Color','w');
%         
%     %end;
% end;




