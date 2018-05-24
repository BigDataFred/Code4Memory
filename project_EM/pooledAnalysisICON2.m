%%
addpath('/home/rouxf/prj/Bham/code/mcode/utils/LFP/');
addpath('/home/rouxf/prj/Bham/code/mcode/utils/');

p2d = '/media/rouxf/rds-share/iEEG_DATA/MICRO/code4MEM/';
origfile = dir([p2d,'*_spikeANDwvfDataAll.mat']);
files1 = dir([p2d,'*_spikeANDwvfData.mat']);

%%
t = [-2 5]; % time window of interest
t = t.*1e3; % convert s to ms
dt = t(1):250:t(2);% bins used for the isi range

cnt = 0;cnt2 = 0;
allSTA = [];
allSTAbp = [];
STApow  = [];
STApowc = [];
STP     = [];
STPc    =[]; 
instFR = [];
selCode = [];
SFC     = [];
Sxx     = [];
Sxxc     = [];
STAen     = [];
STAcu     =[];
STPen     = [];
STPcu     =[];
glFR = [];
PPC     =[];
PLV = [];
ITC = [];
Stf = [];
Ctf = [];
hPL = [];
PLp = [];
ERP = [];
PHI = {};
uCode = [];

%%
for it = 1:length( files1 )
    fprintf([num2str(it),'/',num2str(length(files1))]);
    
    origDat = load([p2d,origfile(it).name]);
    selIx2 = [find(origDat.uSelIx(:,8));find(origDat.uSelIx(:,9))];
    uCode = [uCode ones(1,length(find(origDat.uSelIx(:,8)))) 2*ones(1,length(find(origDat.uSelIx(:,9))))];
    
    for jt = 1:length(selIx2)
        cnt2 = cnt2+1;
        x = origDat.instFR{selIx2(jt)};
        b = origDat.instFR{selIx2(jt)}(dt >=-1e3 & dt <=0);
        z = (x-median(b))./std(b);
        instFR(cnt2,:) = z;
    end;

    dat1 = load([p2d,files1(it).name]);
    
    u1 = find(dat1.uSelIx(:,8) ~=0 & dat1.uSelIx(:,9) ==0);
    u2 = find(dat1.uSelIx(:,9) ~=0 & dat1.uSelIx(:,8) ==0);
    u3 = find(dat1.uSelIx(:,8) ~=0 & dat1.uSelIx(:,9) ~=0);
    selIx = [u1;u2;u3];%1:size(dat1.uSelIx,1);%
    
    selCode = [selCode [ones(1,length(u1)) 2*ones(1,length(u2)) 3*ones(1,length(u3))]];
    
    if ~isempty(selIx)
        
        sesh = files1(it).name(1:max(regexp(files1(it).name,'\d{2}_'))+1);
        
        files2 = dir([p2d,sesh,'_STAandSTPandSTApow.mat']);
        files3 = dir([p2d,sesh,'_spectrumPOWandSFC.mat']);
        files4 = dir([p2d,sesh,'_STPandSTApowCueANDEncoding.mat']);
        files5 = dir([p2d,sesh,'_PPC.mat']);
        files6 = dir([p2d,sesh,'_PLV.mat']);
        files7 = dir([p2d,sesh,'_ITC.mat']);
        files8 = dir([p2d,sesh,'_timeFrequencyData.mat']);
        
        dat2 = load([p2d,files2.name]);
        dat3 = load([p2d,files3.name]);
        dat4 = load([p2d,files4.name]);
        dat5 = load([p2d,files5.name]);
        dat6 = load([p2d,files6.name]);
        dat7 = load([p2d,files7.name]);
        dat8 = load([p2d,files8.name],'spectralData');
        
        for jt = 1:length(selIx)            
            
            cnt = cnt+1;
            
            PHI{cnt} = dat6.phi{selIx(jt)};
            
            x = squeeze(mean(dat8.spectralData.SgmLFP{selIx(jt)},3));
            Stf(cnt,:,:) = x;
            x = squeeze(mean(dat8.spectralData.Cgrm{selIx(jt)},3));
            Ctf(cnt,:,:) = x;
            
            PPC(cnt,:) = dat5.PPC(selIx(jt),:);
            PLV(cnt,:) = dat6.PL{selIx(jt)};
            hPL(cnt) = dat6.hPL(selIx(jt));
            PLp(cnt,:) = dat6.pvalPL{selIx(jt)};
            
            chL = dat7.BFinf(1,selIx(jt))';
            id = unique(dat7.BFinf(1,:));
            ix = find(strcmp(id,chL));
            
            ITC(cnt,:,:) = squeeze(dat7.ITC(ix,:,:));
            
            [b,a] = butter(4,30./1e3,'low');
            x = squeeze(dat7.ERPimg(:,ix,:));
            for it = 1:size(x,1)
                x(it,:) = filtfilt(b,a,x(it,:));
            end;
            
            ERP(cnt,:) = mean(x,1);
            
            x = dat2.allSTA{selIx(jt)};
            y = dat2.allSTAbp{selIx(jt)};
            parfor kt = 1:size(x,1)
                x(kt,:) = x(kt,:) - mean(x(kt,:));
                y(kt,:) = y(kt,:) - mean(y(kt,:));
            end;
            allSTA(cnt,:) = mean(x,1);
            allSTAbp(cnt,:) = mean(y,1);
            
            fx = dat2.allSTP{selIx(jt)}{1}.fx1';
            Snn1 = mean(dat2.allSTP{selIx(jt)}{1}.S1,2);
            STApow(cnt,:) = Snn1;
            
            [Snn1] = overfRegression(Snn1,fx);
            STApowc(cnt,:) = Snn1;
            
            Snn2 = dat2.STApow{selIx(jt)}{1}.S1;
            STP(cnt,:)    =  Snn2;
            
            [Snn2] = overfRegression(Snn2,fx);
            STPc(cnt,:) = Snn2;
            
            SFC(cnt,:) = dat3.spectrumSFC{selIx(jt)}.SFC1;
            Sxx(cnt,:) = mean(dat3.spectrumPOW{selIx(jt)}.S1,2);
            fx = dat3.spectrumPOW{selIx(jt)}.fx1';
            [Sxxc(cnt,:)] = overfRegression(Sxx(cnt,:)',fx);
            
            fx = dat4.STApow1{selIx(jt)}{1}.fx1';
            
            Snn =  dat4.STApow1{selIx(jt)}{1}.S1;
            [Snnc] = overfRegression(Snn,fx);            
            STAcu(cnt,:)     = Snnc;   
            
            Snn =  dat4.STApow2{selIx(jt)}{1}.S1;
            [Snnc] = overfRegression(Snn,fx);
            STAen(cnt,:)     = Snnc;   
            
            Snn =  mean(dat4.allSTP1{selIx(jt)}{1}.S1,2);
            [Snnc] = overfRegression(Snn,fx);
            STPcu(cnt,:)     = Snnc;
            
            Snn =  mean(dat4.allSTP2{selIx(jt)}{1}.S1,2);
            [Snnc] = overfRegression(Snn,fx);
            STPen(cnt,:)     = Snnc;
            
            glFR(cnt) = mean(dat1.glFR{selIx(jt)});
        end;
    end;
    fprintf('\n');
end;

%%
M = mean(instFR(:,2:end-1),1);
SE = std(instFR(:,2:end-1),0,1)./sqrt(size(instFR,1)-1);
figure;
hold on;
h = [];
h(1) = area([0.5 1.5].*1e3,ones(1,2)*max(M+SE),min(M-SE));
h(2) = area([2.5 3.5].*1e3,ones(1,2)*max(M+SE),min(M-SE));
set(h,'FaceColor',[.75 .75 .75],'EdgeColor',[.75 .75 .75],'FaceAlpha',.3,'EdgeAlpha',.3);
jbfill(dt(2:end-1),M-SE,M+SE,[.9 0 0],[.9 0 0],0,.4);
hold on;
plot(dt(2:end-1),M,'r','LineWidth',3);
plot([dt(2) dt(end-1)],[0 0],'k--');
axis tight;
xlim([-250 dt(end-1)]);
box off;
xlabel('Time rel. to cue onset [s]');
ylabel('Norm. firing rate [\sigma]');
set(gca,'XTick',[0 2 4].*1e3);

%%
fx = dat3.spectrumPOW{1}.fx1;

figure;
subplot(121);
a = gca;
plot(fx,log10(Sxx),'Color',[.75 .75 .75]);
subplot(122);
a = [a gca];
hold on;
plot(fx,Sxxc,'Color',[.75 .75 .75]);
plot(fx,mean(Sxxc,1),'k','LineWidth',3);
axis tight;

fx = dat2.STApow{1}{1}.fx1;

figure;
subplot(221);
a = [a gca];
hold on;
plot(fx,log10(STApow),'Color',[.75 .75 .75]);
plot(fx,mean(log10(STApow)),'k','LineWidth',3);
Y = mean(log10(STApow),1)';
X = [ones(size(Y)) log10(fx')];
b = X\Y;
yp = X*b;
plot(fx,yp,'r','LineWidth',3);
subplot(222);
a = [a gca];
plot(fx,STApowc,'b');
axis tight;

subplot(223);
a = [a gca];
hold on;
plot(fx,log10(STP),'Color',[.75 .75 .75]);
plot(fx,mean(log10(STP)),'k','LineWidth',3);
Y = mean(log10(STP),1)';
X = [ones(size(Y)) log10(fx')];
b = X\Y;
yp = X*b;
plot(fx,yp,'r','LineWidth',3);
subplot(224);
a = [a gca];
plot(fx,STPc,'b');
axis tight;

for it = 1:length(a)
    xlabel(a(it),'Frequency [Hz]');
end;

for it = [3 5]
    ylabel(a(it),'Power [log]');
end;

for it = [4 6]
    ylabel(a(it),'Power 1/f adjusted [a.u.]');
end;
axis(a,'tight');
set(a,'box','on');

%%
tAx = linspace(-.5,.5,999);
figure;
hold on;
plot(tAx,mean(allSTA,1),'k','LineWidth',3);
plot(tAx,mean(allSTAbp,1),'g','LineWidth',3);
xlabel('Time [s]');
ylabel('Amplitude [\muV]');
% n = round(size(allSTA,1)/5);
% figure;
% for it = 1:size(allSTA,1)
%     subplot(n,5,it);
%     hold on;
%     plot(tAx,allSTA(it,:),'k');
%     plot(tAx,allSTAbp(it,:),'g','LineWidth',3);
%     
% end;
%%
% fx = dat2.STApow{1}{1}.fx1;
% pow1 = mean(STPc(:,fx>3 & fx <8),2);
% pow2 = mean(STPc(:,fx>12 & fx <24),2);
% 
% ix1 = find(pow1 >=median(pow1));
% ix2 = find(pow1 <median(pow1));
% ix3 = find(pow2 >=median(pow2));
% ix4 = find(pow2 <median(pow2));
% 
% figure;
% subplot(121);
% a = gca;
% hold on;
% plot(ones(1,length(ix1)),glFR(ix1),'ko','MarkerFaceColor','b');
% plot(2*ones(1,length(ix2)),glFR(ix2),'ko','MarkerFaceColor','r');
% plot([.7 1.3],ones(1,2)*median(glFR(ix1)),'k','LineWidth',3);
% plot([1.7 2.3],ones(1,2)*median(glFR(ix2)),'k','LineWidth',3);
% 
% subplot(122);
% a = [a gca];
% hold on;
% plot(ones(1,length(ix3)),glFR(ix3),'ko','MarkerFaceColor','b');
% plot(2*ones(1,length(ix4)),glFR(ix4),'ko','MarkerFaceColor','r');
% plot([.7 1.3],ones(1,2)*median(glFR(ix3)),'k','LineWidth',3);
% plot([1.7 2.3],ones(1,2)*median(glFR(ix4)),'k','LineWidth',3);
% 
% X = [glFR(ix1) glFR(ix2)]';
% G = {[ones(1,length(ix1)) 2*ones(1,length(ix2))]}';
% %boxplot(X,G,'notch','compact');
% 
% set(a,'Xlim',[0 3]);

%%
fx = dat4.STApow1{1}{1}.fx1;

% figure;
% subplot(221);
% plot(fx,STAcu,'Color',[.75 .75 .75]);
% subplot(222);
% plot(fx,STAen,'Color',[.75 .75 .75]);
% subplot(223);
% plot(fx,STPcu,'Color',[.75 .75 .75]);
% subplot(224);
% plot(fx,STPen,'Color',[.75 .75 .75]);

figure;
subplot(121);
D = STPen-STPcu;
M = mean(D,1)
SE = std(D,0,1)./sqrt(size(D,1)-1);
jbfill(fx,M-SE,M+SE,[0 0 0.9],[0 0 0.9],0,.6);
hold on;
plot(fx,M,'b');
ylabel('\theta-power (Encoding - cue)');
axis tight;
box on;
set(gca,'XTick',[5:10:30]);
xlim([-1 31])

subplot(122);
hold on;
pow2 = mean(STPen(:,fx>=3 & fx <=8),2);
pow1 = mean(STPcu(:,fx>=3 & fx <=8),2);
for it = 1:length(pow1)
    plot([1 2],[pow1(it) pow2(it)],'k');
end;
plot(ones(1,length(pow1)), pow1, 'ko','MarkerFaceColor','r','MarkerSize',8);
plot(2*ones(1,length(pow2)), pow2, 'ko','MarkerFaceColor','b','MarkerSize',8);
xlim([0 3]);
set(gca,'XTick',[1 2]);
set(gca,'XTickLabel',{'Enc' 'Cue'});
ylabel('\theta-power');
axis tight;
box on;xlim([0 3]);
%set(gca,'XTickLabelRotation',-45);

%%
tAx = dat8.spectralData.tx;
fAx = dat8.spectralData.fx;
tAx = tAx-1;
ix1 = find(tAx > 0 & tAx <= 2);
ix2 = find(tAx > 2 & tAx <= 4);

Ztf = zeros(size(Stf));
tp = zeros(size(Stf,1),size(Stf,2));
for it = 1:size(Stf,1)
    M = repmat(mean(Stf(it,tAx<0,:),2),[1 size(Stf,2) 1]);
    SD = repmat(std(Stf(it,tAx<0,:),0,2),[1 size(Stf,2) 1]);
    Ztf(it,:,:) = (Stf(it,:,:) - M)./SD;
    tp(it,:) = squeeze(mean(Ztf(it,:,fAx >=3 & fAx <=9),3));
end;

p = [];h = [];
for it = 1:size(tp,2)
    
   [h(it),p(it)] = ttest(tp(:,it));
    
end;
[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p,0.05);

M = mean(tp,1);
SE = std(tp,0,1)./sqrt(size(tp,1)-1);
y = squeeze(mean(Ztf,1))';
%y = (y-min(min(y)))./(max(max(y))-min(min(y)));

figure;
hold on;
imagesc(tAx,fAx,y);
%plot([0 0],[min(fAx) max(fAx)],'w');
%plot([2 2],[min(fAx) max(fAx)],'w');
axis xy;
xlabel('Time rel. to cue onset [s]');
ylabel('Frequency [Hz]');
ca = caxis;
axis tight;
set(gca,'YTick',[5:10:30]);
set(gca,'XTick',[-1:5]);
set(gca,'XTickLabel',{' ' '0' ' ' '2' ' ' '4'});

figure;
cb = colorbar;
caxis(ca);
zlab = get(cb,'YLabel');
set(zlab,'String','z-score [\sigma]');
axis off;
set(cb,'Ticks',[round(ca(1)*10)/10 floor(ca(2)*10)/10]);
set(cb,'TickLength',get(cb,'TickLength').*2);

figure;
subplot(2,1,1);
jbfill(tAx,M-SE,M+SE,'b','b',0,.4);
hold on;
plot(tAx,mean(tp,1),'b','LineWidth',3);
xlim([tAx(1) tAx(end)]);
axis tight;
xlabel('Time rel. to cue onset [s]');
ylabel('3-8Hz power [\sigma]');
set(gca,'XTick',[-1:5]);
set(gca,'XTickLabel',{' ' '0' ' ' '2' ' ' '4'});
box off;

subplot(2,1,2);
tAx = -1:0.025:5;
M = squeeze(mean(mean(ITC(:,3:6,:),2),1))';
SE = squeeze(std(mean(ITC(:,3:6,:),2),0,1))';
jbfill(tAx,M-SE,M+SE,[.9 0 0],[.9 0 0],0,.5);
hold on;
plot(-1:0.025:5,M);
xlim([tAx(1) tAx(end)]);
axis tight;
xlabel('Time rel. to cue onset [s]');
ylabel('3-8Hz ITC [a.u.]');
set(gca,'XTick',[-1:5]);
set(gca,'XTickLabel',{' ' '0' ' ' '2' ' ' '4'});
box off;

%figure;
%plot(dat7.tAx,mean(ERP,1));

%%
fAx = dat5.freqAx;

foi = 1:30;
pf = [];
for it = 1:size(PPC,1)
    
    ix = find(PPC(it,:) == max(PPC(it,:)));
    pf(it) = fAx(ix);
    
end;

M = mean(PPC,1);
SE = std(PPC,0,1)./sqrt(size(PPC,1)-1);

figure;
jbfill(fAx,M-SE,M+SE,[0 0 .9],[0 0 .9],0,.4);
hold on;
plot(1:30,M,'b','LineWidth',3);
xlabel('Frequency [Hz]');
ylabel('Pairwise phase consistency [a.u.]');
box on;
set(gca,'XTick',[5:10:30]);