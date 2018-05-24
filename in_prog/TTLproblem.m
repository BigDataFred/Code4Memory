%%
addpath('/home/rouxf/prj/Bham/code/mcode/helper/logfile_readers/');

params.p = '/media/rouxf/rds-share/Archive/MICRO/P02/log/EM/';
params.fn ='P02_fVSp_SE01_10-Jul-2016_18_45_5180100_LogFile_EMtask.txt';

[dE,dR,b] = getNewLogDataEM2(params);

%%
dt = [32 64];
selIx = find(ttls ==7);

% selIx2 = [];
% for it = 1:length(selIx)
%     if ttls(selIx(it)+1) == 0
%         selIx2 = [selIx2 selIx(it)+1];
%     elseif ttls(selIx(it)+2) ==0
%         selIx2 = [selIx2 selIx(it)+2];        
%     end;
% end;

blockSz = b;%diff(logDat.LogDat1.idx,[],2)+1;%    b;

ix = 1:blockSz(1);
dum = {};
for it = 1:length(blockSz)
    dum{it} = [ix ix+length(ix)];
    if it < length(blockSz)
        ix = [dum{it}(end)+1:dum{it}(end)+blockSz(it+1)];
    end;
end;

trlENC2 = {};
for it = 1:length(dum)
    trlENC2{it} = dum{it}(1:blockSz(it));
end;
trlENC2 = [trlENC2{:}];

%if unique(diff([trlENC' trlENC2'],[],2))~=0;error('wrong trial assignment');end;
selIx= selIx(trlENC2);

x = [];
for it =1:length(TTLix)    
    x(it,:) = dataSamples(TTLix(it)-dt(1):TTLix(it)+dt(2));    
end;

[b,a] = butter(4,300/(Fs/2),'high');% apply low-pass for LFP
x2 = [];
for it = 1:size(x,1)
    [x2(it,:)] = filtfilt(b,a,x(it,:)');
end;

x = x-repmat(mean(x,2),[1 size(x,2)]);
x2 = x2-repmat(mean(x2,2),[1 size(x,2)]);


tx = [-dt:dt]./Fs.*1e3;
ix = find(tx >=0.4 & tx <=1.5);

ix2 = [];
for it = 1:size(x,1)
    [~,ix2(it)] = max(x2(it,ix));
end;

%%
for it = 1:length(TTLix)
    dataSamples(TTLix(it)-dt(1):TTLix(it)+dt(2)) = dataSamples(TTLix(it)-dt(1):TTLix(it)+dt(2)) - x2(it,:);
end;

%%
figure;
imagesc([-dt(1):dt(2)]./Fs.*1e3,1:size(x,1),x2);
axis xy;
caxis([-2 2]);
xlabel('Time rel. to TTL-timestamp [ms]');
ylabel('TTL-pulse #');
cb = colorbar;
zlab = get(cb,'YLabel');
set(zlab,'String','LFP-Amplitude [\muV]');

figure;
subplot(211);
hold on;
plot([-dt(1):dt(2)]./Fs.*1e3,x2);
xlabel('Time rel. to TTL-timestamp [ms]');
axis tight;
ylabel('LFP-Amplitude [\muV]');
subplot(212);
plot([-dt(1):dt(2)]./Fs.*1e3,mean(x2,1),'Color',[.9 0 0],'LineWidth',3);
axis tight;
xlabel('Time rel. to TTL-timestamp [ms]');
ylabel('Average LFP-Amplitude [\muV]');
%[ax,h1,h2] = plotyy([-dt:dt]./Fs.*1e3,x,[-dt:dt]./Fs.*1e3,mean(x,1));
%axis(ax,'tight');
%set(h1,'Color',[.75 .75 .75]);
%set(h2,'Color',[.9 0 0],'LineWidth',3);
%set(ax(1),'YColor',[.75 .75 .75]);
%set(ax(2),'YColor',[.9 0 0]);
%xlabel(ax(1),'Time rel. to TTL-timestamp [ms]');
%ylabel(ax(1),'LFP-Amplitude [\muV]');
%ylabel(ax(2),'Average LFP-Amplitude [\muV]');

figure;
for it = 1:size(x,1);
    hold on;
    plot(tx,x2(it,:));
    plot(tx(ix(ix2(it))),x2(it,ix(ix2(it))),'ro');
    pause;
    clf;
end;
