%%
addpath(genpath('/home/rouxf/tbx/chronux_2_11/'));
%addpath(genpath('/home/rouxf/tbx/eeglab14_1_1b/'));

addpath(genpath('/home/rouxf/tbx/WagenaarMBL/'));
addpath(genpath('/home/rouxf/prj/Bham/code/mcode/utils/'));

%%
[pID]     = 'P09';
[expMode] = 'fVSpEM';

[rpath]   = '/media/rouxf/rds-share/iEEG_DATA/MICRO/';

[sesh]    = dir([rpath,pID,'/',expMode,'/']);
sesh = sesh([sesh(:).isdir]);

[chck]    = regexp({sesh(:).name},'\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}');

x = [];
for it = 1:length(chck);
    x(it) = ~isempty(chck{it}');
end;
[sesh] = {sesh(find(x)).name}';

nChan = 64;

%%
movingwin1 = [1 0.01];
T = movingwin1(1);
W = 1;
TW = T*W;
k = 2*TW-1;

params1              = [];
params1.Fs           = 1e3;
params1.pad          = 2;
params1.fpass        = [0.5 30];
params1.tapers       = [TW k];
params1.trialave     = 0;

movingwin2 = [0.25 0.01];
T = movingwin2(1);
W = 10;
TW = T*W;
k = 2*TW-1;

params2              = [];
params2.Fs           = 1e3;
params2.pad          = 2;
params2.fpass        = [20 170];
params2.tapers       = [TW k];
params2.trialave     = 0;

%%
afoi = 30:2:100;
pfoi = 2:1:20;
pbins = -pi:pi/9:pi;

Fs = 1e3;

b1 = cell(1,length(pfoi));
for it = 1:length(pfoi)
    [b1{it}] = fir1( 3*floor(Fs/(pfoi(it)-1)) ,[pfoi(it)-1 pfoi(it)+1]./(Fs/2),'bandpass');
end;

b2 = cell(1,length(afoi));
for it = 1:length(afoi)
    [b2{it}] = fir1( 3*floor(Fs/(afoi(it)-2)) ,[afoi(it)-2 afoi(it)+2]./(Fs/2),'bandpass');
end;

%%
% f1 = 6.34;
% f2 = f1*2;
% t  = 0:1/Fs:12;
% sig1 = sin(2*pi*f1.*t);
% sig2 = sin(2*pi*f2.*t);
% sig3 = sig1+sig2+(2.75*randn(1,length(t)));
%sig3 = sig1+((sig2).*(sig1+1))+(2.75*randn(1,length(t)));

% n = length(sig3);
% filt1 = filtfilt(b1{pfoi==f1},1,[fliplr(sig3) sig3 fliplr(sig3)]);
% filt2 = filtfilt(b2{afoi==f2},1,[fliplr(sig3) sig3 fliplr(sig3)]);
% %filt1 = eegfilt([fliplr(sig3) sig3 fliplr(sig3)],Fs,f1-2,f1*2);
% %filt2 = eegfilt([fliplr(sig3) sig3 fliplr(sig3)],Fs,f2-5,f2+5);
% filt1 = filt1(n+1:2*n);
% filt2 = filt2(n+1:2*n);
% 
% figure;
% subplot(131)
% hold on;
% plot(t,sig3,'r');
% plot(t,sig1);
% plot(t,sig2,'k');
% axis tight;
% xlim([0 .5]);
% subplot(132)
% hold on;
% plot(t,sig1);
% plot(t,filt1,'r');
% axis tight;
% xlim([0 .5]);
% subplot(133)
% hold on;
% plot(t,sig2);
% plot(t,filt2,'r');
% axis tight;
% xlim([0 .5]);

% phi = zeros(length(pfoi),length(sig3));
% for it = 1:length(pfoi)
%     [phi(it,:)] = angle(hilbert(filtfilt(b1{it},1,sig3)));    
% end;
% amp = zeros(length(pfoi),length(sig3));
% for it = 1:length(afoi)
%     [amp(it,:)] = abs(hilbert(filtfilt(b2{it},1,sig3))).^2;    
% end;
% 
% [PAC] = zeros(length(pfoi),length(afoi),length(pbins));
% for it = 1:length(pfoi)
%     p = phi(it,:);
%     [p,sIx] = sort(p);
%     for jt = 1:length(afoi)
%       a = amp(jt,sIx);      
%       m = zeros(1,length(pbins));
%       for kt = 1:length(pbins)-1
%           ix = find(p >= pbins(kt) & p <pbins(kt+1));
%           m(kt) = mean(a(ix));
%       end;
%       m(end) = m(1);
%       PAC(it,jt,:) = m./sum(m);
%     end;
% end;
% H = -squeeze(sum(PAC.*log(PAC),3));
% n = length(pbins);
% MI = (log(n)-H)./log(n);
% 
% figure;
% pcolor(pfoi,afoi,MI');
% shading interp;

%% concatenate LFP data across sessions
dx = cell(1,nChan);
trlLab = cell(1,nChan);
FR = cell(1,nChan);
rTS = cell(1,nChan);

missIdx2 = [];
hitIdx2 = [];
for jt = 1:length(sesh)
    
    %%
    [p2d] = [rpath,pID,'/',expMode,'/',sesh{jt},'/'];
    [LFPfn] = dir([p2d,pID,'_',expMode,'_',sesh{jt},'_lfpDataStimLockedSegmentedAVGdownsampled.mat']);
    [SPKfn] = dir([p2d,pID,'_',expMode,'_',sesh{jt},'_spkDataStimLockedSegmented.mat']);
    
    %%
    load([p2d,LFPfn.name]);
    spkDat = load([p2d,SPKfn.name]);   
        
    %%
    ntrl = size(LFPavg{1},2);
    
    missIdx2 = find(ismember(sort([missIdx;hitIdx]),missIdx));
    hitIdx2 = find(ismember(sort([missIdx;hitIdx]),hitIdx));
    xx = zeros(ntrl,1);
    xx(hitIdx2) = 1;
    
    %%    
    for it =1:length( LFPavg )
        
        fprintf([num2str(it),'/',num2str(length( LFPavg ))]);
        
        %%        
        x = LFPavg{it};
        
        x = x(:,xx==1);
        M = ones(size(x,1),1)*mean(x,1);
        SD = ones(size(x,1),1)*std(x,0,1);
        
        z = (x-M)./SD;
        z = max(abs(z),[],1);
        x = x(:,z < 4);
        M = M(:,z<4);
        x = x-M;
        
        [s1,t1,f1] = mtspecgramc( gradient( x' )', movingwin1 , params1);
        t1 = t1-5;
        [s2,t2,f2] = mtspecgramc( gradient( x' )', movingwin2, params2  );
        t2 = t2-5;
        s1 = s1(t1>=-1 & t1 <=4,:,:); t1 = t1(t1>=-1 & t1 <=4);
        s2 = s2(t2>=-1 & t2<=4,:,:); t2 = t2(t2>=-1 & t2 <=4);       
        figure;imagesc(t1,f1,squeeze(mean(s1,3))');axis xy;
        figure;imagesc(t2,f2,squeeze(mean(s2,3))');axis xy;
    
        %%
        %x = gradient( x' )';
        for zt = 30%25:32
            fprintf([num2str(zt)])
            x = x';
            [s1,t1,f1] = mtspecgramc( gradient( x' )', movingwin1 , params1);
            [s2,t2,f2] = mtspecgramc( gradient( x' )', movingwin2, params2  );
            figure;plot(f1,squeeze(mean(mean(20*log10(s1),3),1)));
            figure;plot(f2,squeeze(mean(mean(20*log10(s2),3),1)));
            fprintf('\n');
        end;
        
        %%
        for zt = 30%25:32
            fprintf([num2str(zt)])           
            samp = 1:size(x,1);%find(trlTime >=-1 & trlTime <=4);
            MI = zeros(length(pfoi),length(afoi),size(x,2));
            for nt = 1:size(x,2)
                n = length(x(:,nt)');
                phi = zeros(length(pfoi),length(samp));
                for kt = 1:length(pfoi)
                    [s] = filtfilt(b1{kt},1,[fliplr(x(:,nt)') x(:,nt)' fliplr(x(:,nt)') ]);
                    s = angle(hilbert(s(n+1:2*n)));
                    phi(kt,:) = s(samp);
                end;
                
                amp = zeros(length(pfoi),length(samp));
                for kt = 1:length(afoi)
                    [s] = filtfilt(b2{kt},1,[fliplr(x(:,nt)') x(:,nt)' fliplr(x(:,nt)') ]);
                    s = abs(hilbert(s(n+1:2*n))).^2;
                    amp(kt,:) = s(samp);
                end;
                
                [PAC] = zeros(length(pfoi),length(afoi),length(pbins));
                for kt = 1:length(pfoi)
                    p = phi(kt,:);
                    [p,sIx] = sort(p);
                    for lt = 1:length(afoi)
                        a = amp(lt,sIx);
                        m = zeros(1,length(pbins));
                        for mt = 1:length(pbins)-1
                            ix = find(p >= pbins(mt) & p <pbins(mt+1));
                            m(mt) = mean(a(ix));
                        end;
                        m(end) = m(1);
                        PAC(kt,lt,:) = m./sum(m);
                    end;
                end;
                H =squeeze(-sum(PAC.*log(PAC),3));
                n = length(pbins);
                mi = (log(n)-H)./log(n);
                MI(:,:,nt) = mi;
            end;
            figure;
            pcolor(pfoi,afoi,squeeze(mean(MI,3))');
            shading interp;
            fprintf('\n');
        end;
        
        %%

               
        [b] = fir1(3*floor(Fs/1),[1 100]./(Fs/2),'bandpass');
        f = filtfilt(b,1,x);
        
        [b] = fir1(3*floor(Fs/9),[9 11]./(Fs/2),'bandpass');
        f2 = filtfilt(b,1,x);
        
       [ix1] = findpeaks(f);
       [ix2] = findpeaks(-f);
       [ix3] = findpeaks(f2);
       [ix4] = findpeaks(-f2);
       
       sel1 = [];
       for kt = 1:length(ix3.loc)
           [~,sel1(kt)] = min(abs(ix3.loc(kt)-ix1.loc));
       end;
       
       sel2 = [];
       for kt = 1:length(ix4.loc)
           [~,sel2(kt)] = min(abs(ix4.loc(kt)-ix2.loc));
       end;
       ix1.loc = ix1.loc(sel1);
       ix2.loc = ix2.loc(sel2);
       
       dx = diff(ix1.loc);
       delIx = [find(dx <floor(1/90.*Fs))+1];
       ix1.loc(delIx) = [];
       ix2.loc(delIx) = [];
       dx = diff(ix1.loc);
       
       figure;
       plot(f);hold on;plot(f2,'r');
       plot(ix1.loc,f2(ix1.loc),'b^');
       plot(ix2.loc,f2(ix2.loc),'r^');
       
       dt = 500;
       eP = zeros(length(ix2.loc),length(-dt:dt));
       c = 0;
       selIx = [];
       for kt = 1:length(ix2.loc)
           if (ix2.loc(kt) -dt>0) && (ix2.loc(kt)+dt<length(f))
               c = c+1;
               selIx(c) = kt;
               eP(kt,:) = f(ix2.loc(kt)-dt:ix2.loc(kt)+dt);               
           end;           
       end;
       eP = eP(selIx,:)';
       M = ones(size(eP,1),1)*mean(eP,1);
       SD = ones(size(eP,1),1)*std(eP,0,1);
       eP = (eP-M)./SD;
       
       [Sx,fx] = mtspectrumc( eP, params1);
       p = mean(Sx(fx>7 & fx <11,:),1);
       sel1 = find(p>median(p));
       sel2 = find(p<median(p));
       [Sx,tx,fx] = mtspecgramc( gradient(eP')', movingwin2, params2);
       
        %%        
        fprintf('\n');
    end;
    
%     %%
%     BFlab = [];
%     for it = 1:length(chanLab)
%         BFlab = [BFlab;it*ones(8,1)];
%     end;
% 
%     %%
%     for it = 1:length(spkDat.sortedSpikes)
%         spk = spkDat.sortedSpikes{it};
%         trlID = 1:ntrl;
%         trlID = trlID( selIx{BFlab(it)} );
%         dt = -1009:10:4010;
%         n = zeros(length(trlID),length(dt));
%         for kt = 1:length(trlID)
%             spkTs = spk.SpikeTimesSeg(spk.trl == trlID(kt)).*1e3;
%             n(kt,:) = hist(spkTs,dt);
%             rTS{it}{jt}{kt} = spkTs;
%         end;
%         FR{it} = [FR{it};n];        
%     end;
    
end;
nChan = length( chanLab);

dx = dx(1:nChan);
trlLab = trlLab(1:nChan);

%%
n = [];
for it = 1:length( dx )
    [n(it)] = size(dx{it},1);
end;
dx = dx(n>20);
trlLab = trlLab(n>20);
chanLab = chanLab(n>20);
FR = FR(ismember(BFlab,find(n>20)));
rTS = rTS(ismember(BFlab,find(n>20)));
BFlab = BFlab(ismember(BFlab,find(n>20)));
n = n(n>20);

nChan = length( chanLab);



%%
S1 = cell(1,nChan);
S2 = cell(1,nChan);

for it = 1:length( dx )    
    fprintf([num2str(it),'/',num2str(length(dx))]);
    [s1,t1,f1] = mtspecgramc( dx{it}', movingwin1 , params1);
    t1 = t1-5;
    [s2,t2,f2] = mtspecgramc( dx{it}', movingwin2, params2  );
    t2 = t2-5;
    s1 = s1(t1>=-1 & t1 <=4,:,:); t1 = t1(t1>=-1 & t1 <=4);
    s2 = s2(t2>=-1 & t2<=4,:,:); t2 = t2(t2>=-1 & t2 <=4);
    S1{it} = s1;%
    S2{it} = s2;%
    fprintf('\n');
end;
%%
tval1 = [];
tval2 = [];
for it = 1:length( S1 )

    c = squeeze(mean(S1{it}(t1>=0 & t1 <2,:,trlLab{it}==1),1));
    e = squeeze(mean(S1{it}(t1>=2 & t1 <4,:,trlLab{it}==1),1));
    
    for jt = 1:size(c,1)
        [h,p,ci,stats] = ttest(c(jt,:)',e(jt,:)');
        [ tval1(it,jt) ] = stats.tstat;
    end;
    
    c = squeeze(mean(S2{it}(t2>=0 & t2 <2,:,trlLab{it}==1),1));
    e = squeeze(mean(S2{it}(t2>=2 & t2 <4,:,trlLab{it}==1),1));
    for jt = 1:size(c,1)
        [h,p,ci,stats] = ttest(c(jt,:)',e(jt,:)');
        [ tval2(it,jt) ] = stats.tstat;
    end;
    
end;

%%
S1{it} = squeeze(mean(s1(:,:,trlLab{it}==1),3));
    
S2{it} = squeeze(mean(s2(:,:,trlLab{it}==1),3));    
    
%%
[hemLab] = cell( 1 , length( chanLab ) );
for it = 1:length( chanLab )
    hemLab(it) = { chanLab{it}(end) };
end;
[hemLab,sIx] = sort(hemLab);
[ixL] = find(strcmp(hemLab,'L'));
[ixR] = find(strcmp(hemLab,'R'));

AVG1 = zeros(length(t1),length(f1));
for it = 1:length(ixL)
    AVG1 = AVG1 + S1{sIx(ixL(it))};
end;
AVG1 = AVG1./it;

AVG2 = zeros(length(t1),length(f1));
for it = 1:length(ixR)
    AVG2 = AVG2 + S1{sIx(ixR(it))};
end;
AVG2 = AVG2./it;

AVG3 = zeros(length(t2),length(f2));
for it = 1:length(ixL)
    AVG3 = AVG3 + S2{sIx(ixL(it))};
end;
AVG3 = AVG3./it;

AVG4 = zeros(length(t2),length(f2));
for it = 1:length(ixR)
    AVG4 = AVG4 + S2{sIx(ixR(it))};
end;
AVG4 = AVG4./it;

%%
for it = 1:length( chanLab)
    
    Y1 = squeeze(mean( S2{it} ,3));%./squeeze(mean(S2{it}(:,:,trlLab{1}==0),3)));
    Y2 = squeeze(mean( S1{it} ,3));%./squeeze(mean(S1{it}(:,:,trlLab{1}==0),3)));
    
%     m = ones(size(Y1,1),1)*mean(Y1(t1<=0,:),1);
%     Y1 = (Y1-m)./m;
%     
%     m = ones(size(Y1,1),1)*mean(Y2(t2<=0,:),1);
%     Y2 = (Y2-m)./m;
    
    figure;
    subplot(10,1,1:6);
    hold on;
    imagesc(t2,f2,Y1');axis xy;axis tight;xlim([-1 4]);
    plot([0 0],[min(f2) max(f2)],'w');
    plot([2 2],[min(f2) max(f2)],'w');
    title([chanLab{it},'(n= ',num2str(n(it)),')']);
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    subplot(10,1,8:10);
    hold on;
    imagesc(t1,f1,Y2');axis xy;axis tight;xlim([-1 4]);
    plot([0 0],[min(f1) max(f1)],'w');
    plot([2 2],[min(f1) max(f1)],'w');
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
end;

%%
Y1 = AVG1;
Y2 = AVG3;
%m = ones(size(Y1,1),1)*mean(Y1(t1<=0,:),1);
%Y1 = (Y1-m);
  
%m = ones(size(Y2,1),1)*mean(Y2(t2<=0,:),1);
%Y2 = (Y2-m);

figure;
subplot(10,1,1:6);
hold on;
imagesc(t2,f2,Y2');axis xy;axis tight;xlim([-1 4]);
plot([0 0],[min(f2) max(f2)],'w');
plot([2 2],[min(f2) max(f2)],'w');
title(['Left Hemisphere','(n= ',num2str(sum(n(sIx(ixL)))),')']);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
subplot(10,1,8:10);
hold on;
imagesc(t1,f1,Y1');axis xy;axis tight;xlim([-1 4]);
plot([0 0],[min(f1) max(f1)],'w');
plot([2 2],[min(f1) max(f1)],'w');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
    
%%
Y1 = AVG2;
Y2 = AVG4;

figure;
subplot(10,1,1:6);
hold on;
imagesc(t2,f2,Y2');axis xy;axis tight;xlim([-1 4]);
plot([0 0],[min(f2) max(f2)],'w');
plot([2 2],[min(f2) max(f2)],'w');
title(['Right Hemisphere','(n= ',num2str(sum(n(sIx(ixR)))),')']);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
subplot(10,1,8:10);
hold on;
imagesc(t1,f1,Y1');axis xy;axis tight;xlim([-1 4]);
plot([0 0],[min(f1) max(f1)],'w');
plot([2 2],[min(f1) max(f1)],'w');
xlabel('Time (s)');
ylabel('Frequency (Hz)');

%%
pow = cell( 1 , length(S1) );
for it = 1:length( S1 )
    
    pow{it}(:,1) = mean(S1{it}(:,f1 >= 4 & f1 <=12),2);
    pow{it}(:,2) = mean(S1{it}(:,f1 >= 12 & f1 <=30),2);
    pow{it}(:,3) = mean(S2{it}(:,f2 >= 30 & f2 <=90),2);
end;

%%
gw = hanning(10);
fr = zeros(length(FR),length(t1));
for it = 1:length(FR)
    
    x = FR{it};
    x(:,[1 end]) = [];
    x = sum(x,1)./size(x,1)./0.01;
    x = conv(x,gw,'same')./sum(gw);
    fr(it,:) = x;
end;

%%
rXY = cell(1,length(pow));
for it = 1:length( pow )
    
    ix = find(BFlab == it);
    
    x = fr(ix,:)';
    y = pow{it};
    
    if ~isempty(x)
        c = corr(x,y,'Type','Spearman');
    else
        c = NaN(8,3);
    end;
    rXY{it} = c;
    
end;

%%
figure;
hold on;
x = rTS{20};
c = 0;
for it = 1:length(x)
    
    x2 = x{it};
    for jt = 1:length(x2)
        c = c+1;
        s = x2{jt};
        s = s(s>=-1000 & s<=4000);
        
        t = c*ones(1,length(s));
        
        s = [s;s];
        t = [t-.5;t+.5];
        
        line(s,t,'Color','k');
        
    end;
    
end;

%%
s = rTS{11}{1};

n = [];
dt = 3001:4000;
for it = 1:length(s)
    x = s{it};
    x = x(x>=3001 & x <= 4000);
    [n(it,:)] = hist(x,dt);
end;
n = n';

n2 = [];
dt = -999:0;
for it = 1:length(s)
    x = s{it};
    x = x(x>=-999 & x <= 0);
    [n2(it,:)] = hist(x,dt);
end;
n2 = n2';

T = 1;
W = 8;
TW = T*W;
k = 2*TW-1;

params                  = [];
params.pad              = 2;
params.Fs               = 1e3;
params.tapers           = [TW k];
params.fpass            = [25 170];
params.trialave         = 1;

[Sx,fx] = mtspectrumpb( n , params );
[Sx2,fx2] = mtspectrumpb( n2 , params );

figure;plotyy(f2,mean(S2{2}(t2>=3,:),1),f2,mean(S2{2}(t2<=0,:),1));
figure;plotyy(fx,Sx,fx2,Sx2);















