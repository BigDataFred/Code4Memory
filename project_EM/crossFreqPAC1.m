function crossFreqPAC1(pID,expMode)

%%
if nargin == 0
    pID = 'P09';%dataset label
    expMode = 'fVSpEM';
end;

%%
addpath(genpath('/home/rouxf/tbx/chronux_2_11/'));
addpath(genpath('/home/rouxf/tbx/eeglab14_1_1b/functions/sigprocfunc/'));
addpath(genpath('/home/rouxf/prj/Bham/code/mcode/utils/'));

%% parameters for spectral analysis
Fs = 1e3;

T = 391/Fs;% length of time window in s
W = 1/T;% smoothing (+/- 1Hz)
TW = T*W;% time-bandwidth product
k = floor(2*(T*W)-1);% number of tapers
params1              = [];
params1.Fs           = Fs;
params1.pad          = 6;
params1.fpass        = [0.5 30];
params1.tapers       = [TW k];
params1.trialave     = 1;

params              = [];
params.Fs           = Fs;
params.pad          = 4;
params.fpass        = [20 220];
params.trialave     = 1;

movingwin3          = [.05 .001];
T = movingwin3(1);
W = 20;%40;
params5             = params;
params5.tapers      = [T*W 2*(T*W)-1];

movingwin4          = [.1 .001];
T = movingwin4(1);
W = 20;
params6             = params;
params6.tapers      = [T*W 2*(T*W)-1];

%% pre-calculate filter coefficients using finite impulse response filter (linear phase)
bpf = [1 28];
[b1] = fir1(3*floor(Fs/bpf(1)),[ bpf ]./(Fs/2),'bandpass');

bpf = [1 60];
[b2] = fir1(3*floor(Fs/bpf(1)),[ bpf ]./(Fs/2),'bandpass');

bpf = [30 90];
[b3] = fir1(3*floor(Fs/bpf(1)),[ bpf ]./(Fs/2),'bandpass');

bpf = [31];
[b4] = fir1(3*floor(Fs/bpf(1)),[ bpf ]./(Fs/2),'high');

bpf = [30];
[b5] = fir1(3*floor(Fs/bpf(1)),[ bpf ]./(Fs/2),'low');

pfoi = 4:2:20;
fbp = cell(1,length(pfoi));
parfor it = 1:length( pfoi )
    bpf = [pfoi(it)-1 pfoi(it)+1];
    fbp{it} = fir1(3*floor(Fs/bpf(1)),[ bpf ]./(Fs/2),'bandpass');
end;

afoi = 40:4:200;
fba = cell(1,length(afoi));
parfor it = 1:length( afoi )
    bpf = [afoi(it)-2 afoi(it)+2];
    fba{it} = fir1(3*floor(Fs/bpf(1)),[ bpf ]./(Fs/2),'bandpass');
end;

%%
[rpath] = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',pID,'/',expMode,'/'];

%[savePath] = [rpath(1:regexp(rpath,expMode)-1),'res/'];
[savePath] = '/home/rouxf/tmp/figuresLNM/';

chck = dir(savePath);
if isempty(chck)
    mkdir(savePath);
end;

[sesh] = dir(rpath);
[sesh] = sesh([sesh.isdir]==1);
chck = regexp({sesh.name},'\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}');
cnt = 0;
sel = [];
for it = 1:length( chck )
    if ~isempty( chck{it} )
        cnt = cnt+1;
        sel(cnt) = it;
    end;
end;
[sesh] = {sesh(sel).name}'

[p2d] = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',pID,'/',expMode,'/'];

%% load spk and LFP data
for seshIt = 1%:length( sesh)
    
    lfpDat  = [];
    lfpDat = load([p2d,sesh{seshIt},'/',pID,'_',expMode,'_',sesh{seshIt},'_lfpDataStimLockedSegmenteddownsampled.mat']);
        
    %% time range of interest and fS
    lfpDat.trlTime = lfpDat.dsTrlTime;
    [tIx] = find(lfpDat.trlTime >=-1 & lfpDat.trlTime <=4);
    [Fs] = lfpDat.dsFs;
    selTrl = {};
    for it = 1:length(lfpDat.LFPseg)
        selTrl(it) = {find(ismember(lfpDat.trlSel,lfpDat.trlENC))};
    end;
        
    %%
    trlChck = zeros(1,length( selTrl) );
    for it = 1:length( selTrl )
        trlChck(it) = length(selTrl{it});
    end;
    [alphaTrsh] = 0.01/sum(trlChck>=20);
    
    %%
    [PACsig] = zeros(length(lfpDat.chanLab),1);
    pbins = 0:9:360;
    dt = 800;
    dt2 = 195;
    
    %%
    tic;
    for it = 26%1:length(lfpDat.chanLab) % loop over the BFs
        
        fprintf([num2str(it),'/',num2str(length(lfpDat.chanLab)),'\n']);
        
        %%
        %if (length(selTrl{it}) >=20)
            
            %load([p2d,pID,'_',expMode,'_',sesh{seshIt},'_',lfpDat.chanLab{it},'_spectralData.mat'],'pvalPowL','pfL');
            
            %%
            %if (pvalPowL<alphaTrsh) && ( (pfL >= 5) && (pfL <=11) )
                
                %%
                [nsmp,ntrl] = size(lfpDat.LFPseg{it});% get dims of input
                
                %%
                if isempty(gcp('nocreate'))
                    parpool(36,'SpmdEnabled', false)
                end;
                
                %%
                %[thetaTmstmps] = detectThetaBouts(lfpDat.LFPseg{it},Fs,[5 11]);
                thetaTmstmps = cell(1,ntrl);
                for jt = 1:ntrl
                   thetaTmstmps{jt} = [min(lfpDat.trlTime(tIx)) max(lfpDat.trlTime(tIx))];
                end;
                
                thetaIx = cell(1,ntrl);
                for jt = 1:ntrl
                    concat = [];
                    for kt = 1:size(thetaTmstmps{jt},1)
                        x = find( lfpDat.trlTime(tIx) >= thetaTmstmps{jt}(kt,1) &  lfpDat.trlTime(tIx) <= thetaTmstmps{jt}(kt,2) );
                        concat = [concat x];
                    end;
                    thetaIx{jt} = concat;
                end;
                clear concat;
                
                %%
                chck = [];for jt = 1:length(thetaIx);chck(jt) = ~isempty(thetaIx{jt});end;
                thetaIx2 = thetaIx(chck==1);
                
                pbins2 = 0:20:360;
                %phi = zeros(ntrl*length(tIx),length(pfoi));
                %amp = zeros(ntrl*length(tIx),length(afoi));
                phi = zeros(length([thetaIx2{:}]),length(pfoi));
                amp = zeros(length([thetaIx2{:}]),length(afoi));
                for jt = 1:length( thetaIx2 )
                    
                    fprintf([num2str(jt),'/',num2str(length( thetaIx2 ))]);
                    
                    if ~isempty(thetaIx2{jt})
                        x = lfpDat.LFPseg{it}(:,jt)';
                        %ix = (jt-1)*length(tIx)+1:(jt-1)*length(tIx)+length(tIx);
                        if jt == 1
                            ix = 1:length(thetaIx2{jt});
                        else
                            ix = ix(end)+1:ix(end)+length(thetaIx2{jt});
                        end;
                        
                        %d = zeros(length(tIx),length(pfoi));
                        d = zeros(length(thetaIx2{jt}),length(pfoi));
                        tIx2 = tIx(thetaIx2{jt});
                        parfor kt = 1:length(pfoi)
                            ft = filtfilt(fbp{kt},1,[fliplr(x) x fliplr(x)]);
                            ft = hilbert(ft);
                            ft = ft(nsmp+1:2*nsmp);
                            ft = ft(tIx2);
                            d(:,kt) = (angle(ft)+pi).*(180/pi);
                        end;
                        phi(ix,:) = d;
                        
                        %d = zeros(length(tIx),length(afoi));
                        d = zeros(length(thetaIx2{jt}),length(pfoi));
                        tIx2 = tIx(thetaIx2{jt});
                        parfor kt = 1:length(afoi)
                            ft = filtfilt(fba{kt},1,[fliplr(x) x fliplr(x)]);
                            ft = hilbert(ft);
                            ft = ft(nsmp+1:2*nsmp);
                            ft = ft(tIx2);
                            d(:,kt) = abs(ft).^2;
                        end;
                        amp(ix,:) = d;
                    end;
                    fprintf('\n');
                end;
                clear d;
                
                X = zeros(length(afoi),length(pfoi),length(pbins2));
                for kt = 1:length( pfoi )
                    p = phi(:,kt);
                    for lt = 1:length( afoi )
                        a = amp(:,lt);
                        parfor mt = 1:length(pbins2)-1
                            X(lt,kt,mt) = mean(a(p>=pbins2(mt) & p < pbins2(mt+1)));
                        end;
                        X(lt,kt,end) = X(lt,kt,1);
                        X(lt,kt,:) = X(lt,kt,:)./sum(X(lt,kt,:));
                    end;
                end;
                H = -squeeze(sum(X.*log(X),3));clear X;
                n = length(pbins2);
                MI = (log(n)-H)./log(n);clear H;
                
                chck = [];for jt = 1:length(thetaIx);chck(jt) = ~isempty(thetaIx{jt});end;
                thetaIx2 = thetaIx(chck==1);
                nrand = 200;
                MI2 = zeros(nrand,length(afoi),length(pfoi));
                parfor nt = 1:nrand
                    fprintf([num2str(nt),'/',num2str(nrand)]);
                    
                    trlp = 1:length( thetaIx2 );
                    %ixa = 1:length(tIx);
                    ixa = [];
                    
                    %shuffle the amplitude time series against the phase time series
                    %amp2 = zeros(ntrl*length(tIx),length(afoi));% not this could also be phase instead
                    %phi2 = zeros(ntrl*length(tIx),length(afoi));% not this could also be phase instead
                    amp2 = zeros(length([thetaIx2{:}]),length(afoi));% not this could also be phase instead
                    phi2 = zeros(length([thetaIx2{:}]),length(pfoi));% not this could also be phase instead
                    
                    f = 0;
                    while f<1
                        rIx = trlp(randperm(length(trlp)))';% pick random order of trials
                        if (sum(diff([trlp' rIx],[],2)==0)==0)
                            f = 1;
                        end;
                    end;
                    
                    ixa = 1:length(thetaIx2{rIx(1)});
                    for jt = 1:length( thetaIx2 )
                        
                        if ~isempty(thetaIx2{rIx(jt)})
                            %ix = (rIx(jt)-1)*length(tIx)+1:(rIx(jt)-1)*length(tIx)+length(tIx);
                            ix = length([thetaIx2{1:rIx(jt)-1}])+1:length([thetaIx2{1:rIx(jt)-1}])+length([thetaIx2{rIx(jt)}]);
                            
                            for kt = 1:length(afoi)
                                amp2(ixa,kt) = amp(ix,kt);% amp
                            end;
                            for kt = 1:length(pfoi)
                                phi2(ixa,kt) = phi(ix,kt);% phase
                            end;
                            %ixa = ixa+length(tIx);
                            if jt < length( thetaIx2 )
                                ixa = ixa(end)+1:ixa(end)+length([thetaIx2{rIx(jt+1)}]);
                            end;
                        end;
                    end;
                    
                    X = zeros(length(afoi),length(pfoi),length(pbins2));
                    for kt = 1:length( pfoi )
                        cf = 1;%randperm(2);
                        if cf(1) ==1
                            p = phi(:,kt);%orig phase
                        else
                            p = phi2(:,kt);% random phase
                        end;
                        for lt = 1:length( afoi )
                            if cf(1)==1
                                a = amp2(:,lt);% random amp
                            else
                                a = amp(:,lt);% orig amp
                            end;
                            for mt = 1:length(pbins2)-1
                                X(lt,kt,mt) = mean(a(p>=pbins2(mt) & p < pbins2(mt+1)));
                            end;
                            X(lt,kt,end) = X(lt,kt,1);
                            X(lt,kt,:) = X(lt,kt,:)./sum(X(lt,kt,:));
                        end;
                    end;
                    H = -squeeze(sum(X.*log(X),3));
                    n = length(pbins2);
                    MI2(nt,:,:) = (log(n)-H)./log(n);
                    fprintf('\n');
                end;
                clear amp phi;
                
                pval = zeros(size(MI));
                for jt = 1:size(MI,1)
                    parfor kt = 1:size(MI,2)
                        ref = squeeze( MI2(:,jt,kt) );
                        pval(jt,kt) = length(find(ref >= MI(jt,kt)))/length(ref);
                    end;
                end;
                clear M2;
                
                [trsh] = (pval < alphaTrsh);
                
                %%
                if any( any( trsh ))
                    
                    PACsig(it) = 1;
                    %[pF] = unique([min(pfoi(find(mean(MI.*trsh,1)~=0))) max(pfoi(find(mean(MI.*trsh,1)~=0)))]);
                    %if length(pF)<2
                        %pF = [pF-1 pF+1];
                        pF = [4 11];
                    %end;
                    %[aF] = unique([min(afoi(find(mean(MI.*trsh,2)~=0))) max(afoi(find(mean(MI.*trsh,2)~=0)))]);
                    %if length(aF)<2
                        %aF = [aF-10 aF+10];
                        aF = [60 80];
                    %end;
                    
%                     saveName = [pID,'_',expMode,'_',sesh{seshIt},'_',lfpDat.chanLab{it},'_comodulogramDat.mat'];
%                     save([savePath,saveName],'MI','afoi','pfoi','pval','alphaTrsh');
                    
                    %%
                    h = [];
                    h(1) = figure;
                    subplot(5,1,1:4);
                    hold on;
                    pcolor(pfoi,afoi,MI);
                    cb = colorbar;
                    axis tight;
                    shading interp;lighting phong;
                    xlabel('Frequnecy for Phase (Hz)');
                    ylabel('Frequency for Amplitude (Hz)');
                    zlab = get(cb,'YLabel');
                    set(zlab,'String','Modulation Index');
                    caxis([0 5e-3]);
                    title([lfpDat.chanLab{it},': Peak-Frequency ',num2str(pF(1)),':',num2str(pF(2)),'Hz']);
                    subplot(5,1,5);
                    plot(pfoi,mean(MI,1),'b','LineWidth',3);xlim([pfoi(1) pfoi(end)]);
                    xlabel('Frequnecy for Phase (Hz)');
                    ylabel('MI');
                    
                    h(2) = figure;
                    subplot(5,1,1:4);
                    hold on;
                    pcolor(pfoi,afoi,MI.*trsh);
                    cb = colorbar;
                    axis tight;
                    shading interp;lighting phong;
                    xlabel('Frequency (Hz)');
                    ylabel('Frequency (Hz)');
                    zlab = get(cb,'YLabel');
                    set(zlab,'String','Modulation Index');
                    caxis([0 5e-3]);
                    title([lfpDat.chanLab{it},': Peak-Frequency ',num2str(pF(1)),':',num2str(pF(2)),'Hz']);
                    subplot(5,1,5);
                    plot(pfoi,mean(MI.*trsh,1),'b','LineWidth',3);xlim([pfoi(1) pfoi(end)]);
                    xlabel('Frequnecy for Phase (Hz)');
                    ylabel('MI');
                    
                    for kt = 1:length(h)
                        saveName = [pID,'_',expMode,'_',sesh{seshIt},'_',lfpDat.chanLab{it},'_comodulogramDat',num2str(kt),'.fig'];
                        saveas(h(kt),[savePath,saveName],'fig');
                    end;
                    clear MI;
                    
                    %%
                    if isempty(gcp('nocreate'))
                        parpool(36,'SpmdEnabled', false)
                    end;
                    
                    %% filter coeffs
                    bpf = [ pF ];
                    [b6] = fir1(3*floor(Fs/bpf(1)),[ bpf ]./(Fs/2),'bandpass');
                    
                    bpf = [ aF ];
                    [b7] = fir1(3*floor(Fs/bpf(1)),[ bpf ]./(Fs/2),'bandpass');
                    
                    %% gamma event trough detection & LFP alignment
                    gP = zeros(ntrl,length(tIx));
                    gE = zeros(ntrl,length(tIx));
                    parfor jt = 1:ntrl
                        [x] = lfpDat.LFPseg{it}(:,jt)';
                        x  = filtfilt(b7,1,[fliplr(x) x fliplr(x)]);% filter LFP to obtain gamma-signal
                        x2 = abs(hilbert(x)).^2;
                        x = x(nsmp+1:nsmp*2);
                        x2 = x2(nsmp+1:nsmp*2);
                        x = x(tIx);
                        x2 = x2(tIx);
                        gP(jt,:) = x;
                        gE(jt,:) = x2;
                    end;
                    M = mean(gE(:));% session mean of gamma activity
                    SD = std(gE(:));% session SD of gamma activity
                    trsh = M+2*SD;
                    
                    ptA1 = zeros(35e3,length(-dt2:dt2));
                    ptA2 = zeros(35e3,length(-dt2:dt2));
                    ptA3 = zeros(35e3,length(-dt2:dt2));
                    cnt2 = 0;
                    phi = [];
                    for jt = 1:ntrl
                        
                        [selIx] = findpeaks( -gP(jt,:) );% detect gamma troughs
                        selIx = selIx.loc;
                        selIx([1 end]) = [];
                        selIx = selIx( gE(jt,selIx)>trsh);%apply treshold
                        
                        y1 = lfpDat.LFPseg{it}(:,jt)';% raw LFP
                        y2 = filtfilt(b6,1,[fliplr(y1) y1 fliplr(y1)]);% filter LFP to obtain theta-signal
                        y3 = filtfilt(b7,1,[fliplr(y1) y1 fliplr(y1)]);% filter LFP to obtain gamma-signal
                        y3 = abs(hilbert(y3)).^2;
                        y4 = angle(hilbert(y2));
                        
                        y2 = y2(nsmp+1:nsmp*2);
                        y3 = y3(nsmp+1:nsmp*2);
                        y4 = y4(nsmp+1:nsmp*2);
                        
                        y1 = y1(tIx);
                        y2 = y2(tIx);
                        y3 = y3(tIx);
                        y4 = y4(tIx);
                        
                        phi = [phi y4(selIx) ];
                        
                        % gamma trough triggered LFP
                        for kt = 1:length(selIx)
                            if (selIx(kt)-dt2 > 0) && (selIx(kt)+dt2<length(tIx))
                                cnt2 = cnt2+1;
                                ptA1(cnt2,:) = y1(selIx(kt)-dt2:selIx(kt)+dt2);%raw LFP
                                ptA2(cnt2,:) = y2(selIx(kt)-dt2:selIx(kt)+dt2);% theta filtered LFP
                                ptA3(cnt2,:) = y3(selIx(kt)-dt2:selIx(kt)+dt2);% gamma envelope of LFP
                            end;
                        end;
                        
                    end;
                    ptA1(cnt2+1:end,:) = [];
                    ptA2(cnt2+1:end,:) = [];
                    ptA3(cnt2+1:end,:) = [];
                    
                    [sIx1] = find(sign(ptA2(:,find(-dt2:dt2==0)))==-1);
                    [sIx2] = find(sign(ptA2(:,find(-dt2:dt2==0)))==1);
                    
                    x1 = ptA3(sIx1,find(-dt2:dt2==0));
                    x1 = (x1-M)./SD;
                    x2 = ptA3(sIx2,find(-dt2:dt2==0));
                    x2 = (x2-M)./SD;
                    
                    ix = find(ptA3(:,find(-dt2:dt2==0)) > trsh);
                    
                    ptA1 = ptA1(ix,:);
                    ptA2 = ptA2(ix,:);
                    ptA3 = ptA3(ix,:);
                    
                    ptA1 = ptA1-(ones(size(ptA1,2),1)*mean(ptA1,2)')';
                    ptA2 = ptA2-(ones(size(ptA2,2),1)*mean(ptA2,2)')';
                    ptA3 = ptA3-(ones(size(ptA3,2),1)*mean(ptA3,2)')';
                    
                    phi = (phi+pi).*(180/pi);
                    
                    [n,x] = hist(phi,0:40:360);
                    
                    [sIx1] = find(sign(ptA2(:,find(-dt2:dt2==0)))==-1);
                    [sIx2] = find(sign(ptA2(:,find(-dt2:dt2==0)))==1);
                    
                    xc1 = zeros(size(ptA2,1),dt2*2+1);
                    xc2 = zeros(size(ptA2,1),dt2*2+1);
                    parfor jt = 1:size(ptA2,1)
                        [xc1(jt,:),~] = xcorr(ptA2(jt,:),dt2,'coeff');
                        [xc2(jt,:),~] = xcorr(ptA2(jt,:),ptA3(jt,:),dt2,'coeff');
                    end;
                    [lag] = -dt2:dt2;
                    
                    M = mean(xc2,1);
                    SE = std(xc2,0,1)/sqrt(size(xc2,1)-1);
                    
                    dum = zeros(size(ptA3));
                    parfor jt = 1:size(ptA3,1)
                        dum(jt,:) = (ptA3(jt,:)-min(ptA3(jt,:)))./(max(ptA3(jt,:))-min(ptA3(jt,:)));
                    end;
                    
                    [Sp,fp] = mtspectrumc( gradient(dum)', params1 );
                    clear dum;
                    
                    readme= {'ptA1: raw LFP','ptA2:theta LFP','ptA3:gamma LFP','xc1:autoCorr theta LFP','xc2:crossCorr thetaLFP-gammaPow','M: mean xc2','SE: standard error xc2'};
                    
                    saveName = [pID,'_',expMode,'_',sesh{seshIt},'_',lfpDat.chanLab{it},'_gammmaTroughTriggeredLFP.mat'];
                    save([savePath,saveName],'ptA1','ptA2','ptA3','dt2','lag','xc1','xc2','x1','x2','M','SE','x','n','fp','Sp','readme');
                    
                    %%
                    h = [];
                    h(1) = figure;
                    hold on;
                    bar([1 2],[mean(x1) mean(x2)]);
                    plot([1 1],[mean(x1) mean(x1)+std(x1)/sqrt(length(x1)-1)],'k');
                    plot([2 2],[mean(x2) mean(x2)+std(x2)/sqrt(length(x2)-1)],'k');
                    plot([.9 1.1],ones(1,2)*mean(x1)+std(x1)/sqrt(length(x1)-1),'k');
                    plot([1.9 2.1],ones(1,2)*mean(x2)+std(x2)/sqrt(length(x2)-1),'k');
                    xlim([0 3]);
                    set(gca,'XTick',1:2);
                    set(gca,'XTickLabel',{[num2str(pF(1)),'-',num2str(pF(2)),'Hz Trough'],[num2str(pF(1)),'-',num2str(pF(2)),'Hz peak']});
                    ylabel([num2str(aF(1)),'-',num2str(aF(2)),'Hz-power [\sigma]']);
                    
                    h(2) = figure;
                    subplot(4,1,1:3);
                    imagesc(-dt2:dt2,1:size(ptA1,1),ptA1([sIx1;sIx2],:));
                    axis xy;
                    xlabel('Time (ms)');
                    ylabel('event#');
                    subplot(4,1,4);
                    hold on;
                    [ax,h1,h2] = plotyy(-dt2:dt2,mean(ptA1,1),-dt2:dt2,mean(ptA2,1));
                    %plot(ax(1),-dt2:dt2,mean(ptA2(sIx1,:),1),'m');
                    %plot(ax(1),-dt2:dt2,mean(ptA2(sIx2,:),1),'g');
                    axis(ax,'tight');
                    axis(ax,'off');
                    
                    h(3) = figure;
                    bar([0:40:760],[n n],'c');
                    axis tight;xlim([0 720]);
                    xlabel([num2str(pF(1)),'-',num2str(pF(2)),'Hz Phase (deg)']);
                    ylabel('Counts');
                    title([num2str(aF(1)),'-',num2str(aF(2)),'Hz peaks']);
                    
                    h(4) = figure;
                    hold on;
                    plot(lag,mean(xc1,1),'k');
                    jbfill(lag,M-SE,M+SE,'r','r',1,.5);
                    hold on;
                    plot(lag,mean(xc2,1),'r');
                    xlabel('Lag (ms)');
                    ylabel('Cross-Correlation Coeff.');
                    
                    h(5) = figure;
                    plot(fp,Sp);axis tight;
                    xlabel('Frequency (Hz)');
                    ylabel('Power');
                    
                    for kt = 1:length(h)
                        saveName = [pID,'_',expMode,'_',sesh{seshIt},'_',lfpDat.chanLab{it},'_gammmaTroughTriggeredLFP',num2str(kt),'.fig'];
                        saveas(h(kt),[savePath,saveName],'fig');
                    end;
                    clear ptA* x1 x2 xc1 xc2 lag M SE readme;
                    
                    %%
                    if isempty(gcp('nocreate'))
                        parpool(36,'SpmdEnabled', false)
                    end;
                    
                    %% theta events trough detection & LFP alignment
                    %avg = zeros(ntrl,205,length(pbins));
                    avg = zeros(ntrl,length(afoi),length(pbins));
                    ptA4 = zeros(20e3,length(-dt:dt));
                    ptA5 = zeros(20e3,length(-dt:dt));
                    cnt = 0;
                    for jt = 1:ntrl
                        fprintf([num2str(jt),'/',num2str(ntrl),'\n']);
                        
                        x   = lfpDat.LFPseg{it}(:,jt)';% raw signal
                        
                        % apply filters
                        ft  = filtfilt(b6,1,[fliplr(x) x fliplr(x)]);% narrow-band theta filtered LFP
                        ft2  = filtfilt(b7,1,[fliplr(x) x fliplr(x)]);% gamma filtered LFP
                        ft3  = filtfilt(b1,1,[fliplr(x) x fliplr(x)]);% wide-band theta filtered LFP
                        
                        % hilbert transform -> analytical phase & amp
                        ft2 = abs(hilbert(ft2)).^2;% power envelope of gamma filtered LFP
                        phi2 = angle(hilbert(ft));% analytical phase of theta filtered LFP
                        
                        %truncate padding
                        ft = ft(nsmp+1:nsmp*2);
                        ft2 = ft2(nsmp+1:nsmp*2);
                        ft3 = ft3(nsmp+1:nsmp*2);
                        phi2 = phi2(nsmp+1:nsmp*2);
                        
                        % extract samples of interest
                        x = x(tIx);
                        ft = ft(tIx);
                        ft2 = ft2(tIx);
                        ft3 = ft3(tIx);
                        phi2 = phi2(tIx);
                        
                        [phWBi,phNBi,ixTWB,ixTNB,ixPhWB,ixPhNB] = interpPhasel(ft3,ft);% use both a wide band and narrow band filtered LFP
                        phWBi = phWBi(ismember(ixPhWB,thetaIx{jt}));
                        phNBi = phNBi(ismember(ixPhNB,thetaIx{jt}));
                        ixPhWB = ixPhWB(ismember(ixPhWB,thetaIx{jt}));
                        ixPhNB = ixPhNB(ismember(ixPhNB,thetaIx{jt}));
                        ixTWB = ixTWB(ismember(ixTWB,thetaIx{jt}));
                        ixTNB = ixTNB(ismember(ixTNB,thetaIx{jt}));
                        
                        %         phi2 = phi2(ix(1,1):ix(end,end));
                        %         phi3 = (phi2+pi).*(180/pi);
                        
                        % theta trough triggered LFP
                        for kt = 1:length(ixTNB)% loop across troughs
                            if (ixTNB(kt)-dt > 0) && (ixTNB(kt)+dt < length(x))
                                cnt = cnt+1;
                                d = x(ixTNB(kt)-dt:ixTNB(kt)+dt);
                                ptA4(cnt,:) = d;% aligh LFP to theta troughs
                                d = ft2(ixTNB(kt)-dt:ixTNB(kt)+dt);
                                ptA5(cnt,:) = d;% aligh gamma power to theta troughs
                            end;
                        end;
                        
                        
                        %                 [S1,t1,f1] = mtspecgramc( gradient(lfpDat.LFPseg{it}(:,jt)')',movingwin3, params5);
                        %                 S1 = S1(t1-5>=-1 & t1-5<=4,:);
                        
                        y  = lfpDat.LFPseg{it}(:,jt);
                        for mt = 1:length(afoi)
                            ft = filtfilt(fba{mt},1,[fliplr(y) y fliplr(y)]);
                            ft = abs(hilbert(ft)).^2;
                            ft = ft(nsmp+1:nsmp*2);
                            ft = ft(tIx);
                            ft = ft(ixPhNB);
                            
                            parfor lt = 1:length(pbins)-1
                                %avg(jt,:,lt) = avg(jt,:,lt)+mean( S1( phNBi>=pbins(lt) & phNBi < pbins(lt+1),: ) ,1);
                                avg(jt,mt,lt) = avg(jt,mt,lt)+mean( ft( phNBi>=pbins(lt) & phNBi < pbins(lt+1)) );
                            end;
                            %avg(jt,:,end) = avg(jt,:,1);
                            avg(jt,mt,end) = avg(jt,mt,1);
                        end;
                        
                        fprintf('\n');
                    end;
                    clear S1 t1;
                    
                    ptA4(cnt+1:end,:) = [];
                    ptA5(cnt+1:end,:) = [];
                    
                    saveName = [pID,'_',expMode,'_',sesh{seshIt},'_',lfpDat.chanLab{it},'_thetaTroughTriggeredLFP.mat'];
                    readme = {'ptA4: raw LFP aligned to theta trough' 'ptA5: gamma Power aligned to theta trough'};
                    save([savePath,saveName],'dt','ptA4','ptA5','readme');
                    %clear ptA5;
                    
                    %%
                    M = repmat(nanmean(nanmean(avg,3),1),[size(avg,1) 1 size(avg,3)]);
                    %x = squeeze(mean(mean(avg(:,f>=50 & f<=100,:),2),3));
                    avg = log(avg)-log(M);
                    
                    %%
                    [S2,t2,f2] = mtspecgramc(gradient(ptA4)',movingwin4,params6);
                    %[S2,t2,f2] = mtspecgramc(gradient(ptA4)',movingwin3,params5);
                    M = ones(size(S2,1),1)*mean(S2,1);
                    S2 = log(S2)-log(M);
                    %     [S1,t,f] = mtspecgramc(gradient(ptA4(x<median(x),:))',movingwin2,params2);
                    %     M = ones(size(S1,1),1)*mean(S1,1);
                    %     S1 = log(S1)-log(M);
                    %     [S2,t,f] = mtspecgramc(gradient(ptA4(x>median(x),:))',movingwin2,params2);
                    %     M = ones(size(S2,1),1)*mean(S2,1);
                    %     S2 = log(S2)-log(M);
                    
                    saveName = [pID,'_',expMode,'_',sesh{seshIt},'_',lfpDat.chanLab{it},'_thetaTroughTriggeredPower.mat'];
                    %save([savePath,saveName],'f1','avg','pbins','S2','t2','f2');
                    save([savePath,saveName],'afoi','avg','pbins','S2','t2','f2');
                    
                    %%
                    h = [];
                    h(1)=figure;
                    hold on;
                    %pcolor(pbins,f1,(squeeze(nanmean(avg,1))));shading interp;
                    pcolor(pbins,afoi,(squeeze(nanmean(avg,1))));shading interp;
                    %plot([180 180],[min(f1) max(f1)],'w--');axis tight;%caxis([-.5 .5]);
                    plot([180 180],[min(afoi) max(afoi)],'w--');axis tight;%caxis([-.5 .5]);
                    set(gca,'XTick',[0:90:360]);
                    xlabel([num2str(pF(1)),'-',num2str(pF(2)),'Hz-phase (deg)']);
                    ylabel('Frequency (Hz)');
                    cb = colorbar;
                    zlab = get(cb,'YLabel');
                    set(zlab,'String',['Spectral Power (dB)']);
                    
                    h(2)=figure;
                    subplot(4,1,1:3);
                    imagesc(t2.*1e3-dt,f2,S2');axis xy;
                    xlim([-300 300]);
                    ylabel('Frequency (Hz)');
                    subplot(4,1,4);
                    [ax,h1,h2] = plotyy(-dt:dt,mean(ptA4,1),t2.*1e3-dt,mean(S2(:,f2>=aF(1) & f2<=aF(2)),2));
                    hold(ax(1),'on');
                    plot(ax(1),-dt:dt,mean(ptA5,1)-mean(mean(ptA5,1)),'k','LineWidth',3);
                    set([h1 h2],'LineWidth',3);
                    box(ax(1),'off');
                    axis(ax,'tight');
                    set(ax,'Xlim',[-300 300]);
                    caxis([-.1 .1]);
                    xlabel('Time (ms)');
                    ylabel(ax(1),'LFP-amplitude (\muV)');
                    ylabel(ax(2),'Power (dB)');
                    
                    for kt = 1:length(h)
                        saveName = [pID,'_',expMode,'_',sesh{seshIt},'_',lfpDat.chanLab{it},'_thetaTroughTriggeredPower',num2str(kt),'.fig'];
                        saveas(h(kt),[savePath,saveName],'fig');
                    end;
                    clear f1 avg t2 f2 S2 ptA4;
                    
                end;
            %end;
        %end;
        close all;
        
    end;
    
    %%
    saveName = [pID,'_',expMode,'_',sesh{seshIt},'_sigPAClist.mat'];
    save([savePath,saveName],'PACsig');
    toc;
    
end;






