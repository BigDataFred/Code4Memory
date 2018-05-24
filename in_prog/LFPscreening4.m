%%
addpath(genpath('/home/rouxf/tbx/chronux_2_11/'));
addpath(genpath('/home/rouxf/prj/Bham/code/mcode/utils/'));

%% parameters for spectral analysis (I)
[Fs] = 1e3;

%lower frequencies
movingwin1 = [1 0.01];
T = movingwin1(1);% length of time window in s
W = 1;% smoothing (+/- 1Hz)
TW = T*W;% time-bandwidth product
k = 2*TW-1;% number of tapers

params1              = [];
params1.Fs           = Fs;
params1.pad          = 2;
params1.fpass        = [0.5 30];
params1.tapers       = [TW k];
params1.trialave     = 1;

%lower frequencies
movingwin2 = [0.25 0.001];
T = movingwin2(1);% length of time window in s
W = 10;% smoothing (+/- 10Hz)
TW = T*W;% time-bandwidth product
k = 2*TW-1;% number of tapers

params2              = [];
params2.Fs           = Fs;
params2.pad          = 2;
params2.fpass        = [20 220];
params2.tapers       = [TW k];
params2.trialave     = 1;

%% pre-calculate filter coefficients using finite impulse response filter (linear phase)
bpf = [4 12];
[b1] = fir1(3*floor(Fs/bpf(1)),[ bpf ]./(Fs/2),'bandpass');

bpf = [1 60];
[b2] = fir1(3*floor(Fs/bpf(1)),[ bpf ]./(Fs/2),'bandpass');

bpf = [60 90];
[b3] = fir1(3*floor(Fs/bpf(1)),[ bpf ]./(Fs/2),'bandpass');

pfoi = 2:2:20;
fbp = cell(1,length(pfoi));
for it = 1:length( pfoi )
    bpf = [pfoi(it)-1 pfoi(it)+1];
    fbp{it} = fir1(3*floor(Fs/bpf(1)),[ bpf ]./(Fs/2),'bandpass');
end;

afoi = 40:4:200;
fba = cell(1,length(afoi));
for it = 1:length( afoi )
    bpf = [afoi(it)-2 afoi(it)+2];
    fba{it} = fir1(3*floor(Fs/bpf(1)),[ bpf ]./(Fs/2),'bandpass');
end;

%%
if isempty(gcp('nocreate'))
    parpool(36,'SpmdEnabled', false)
end;

%%
[rpath] = '/media/rouxf/rds-share/iEEG_DATA/MICRO/';

%% load spk and LFP data
[p2d]       = '/media/rouxf/rds-share/iEEG_DATA/MICRO/P09/fVSpEM/2017-08-28_14-20-18/';
[spkDat]    = load([p2d,'P09_fVSpEM_2017-08-28_14-20-18_spkDataStimLockedSegmented.mat']);
[lfpDat]    = load([p2d,'P09_fVSpEM_2017-08-28_14-20-18_lfpDataStimLockedSegmenteddownsampled.mat']);

%% time range of interest and fS
[tIx] = find(lfpDat.trlTime >=-1 & lfpDat.trlTime <=4);
[Fs] = lfpDat.Fs;

%% parameters for spectral analysis (II)
T = floor(length(tIx)/Fs);
W = 1;%smoothing (+/- 1Hz)
TW = T*W;
k = 2*TW-1;

params3              = [];
params3.Fs           = Fs;
params3.pad          = 2;
params3.fpass        = [0.5 30];
params3.tapers       = [TW k];
params3.trialave     = 1;
params3.err          = [2 0.001];%jacknife error

T = floor(length(tIx)/Fs);
W = 4;%smoothing (+/- 4Hz)
TW = T*W;
k = 2*TW-1;

params4              = [];
params4.Fs           = Fs;
params4.pad          = 2;
params4.fpass        = [20 100];
params4.tapers       = [TW k];
params4.trialave     = 1;
params4.err          = [2 0.001];%jacknife error

%% extract BF labels and ids
ix = regexp(lfpDat.chanLab,'\d{1}');
ix = [ix{:}]-1;
[BFlab] = cell(1,length(ix));
for it = 1:length(ix)
    BFlab(it) = { lfpDat.chanLab{it}(1:ix(it)) };
end;
[BFid] = unique(BFlab)';

%% estimate median LFP
[AVG] = cell(1,length(BFid));
for it = 1:length(BFid)
    
    fprintf([num2str(it),'/',num2str(length(BFid))]);
    
    %% trial indexes: hits vs misses
    ntrl = size(AVG,2);
    
    missIdx2 = find(ismember(sort([lfpDat.missIdx;lfpDat.hitIdx]),lfpDat.missIdx));
    hitIdx2 = find(ismember(sort([lfpDat.missIdx;lfpDat.hitIdx]),lfpDat.hitIdx));
    xx = zeros(ntrl,1);
    xx(hitIdx2) = 1;% keep only hits
    
    %%
    ix = find( strcmp(BFlab,BFid{it}) );
    [nsmp,ntrl] = size( lfpDat.LFPseg{ix(1)}(:,xx==1));
    
    %     %% eigenvalue decomposition
    %     X = zeros(length(ix),nsmp,ntrl);
    %     for jt = 1:length(ix)
    %         X(jt,:,:) = lfpDat.LFPseg{ix(jt)}(:,xx==1);
    %     end;
    %
    %     AVG{it} = zeros(nsmp,ntrl);
    %     for jt = 1:ntrl
    %
    %         x = squeeze(X(:,:,jt));
    %         m = mean(x,2)*ones(1,size(x,2));
    %         x = x-m;
    %         x = x';
    %         R = x'*x;
    %
    %         [v,d] = eig(R);
    %         [~,sIx] = sort(diag(d));
    %         v = v(sIx,:);
    %         sIx = flipud(sIx);
    %         v = v(:,sIx);
    %         y = (v'*x')';
    %         y = y(:,end)';
    %         AVG{it}(:,jt) = y;
    %     end;
    
    %% compute the correlation between microwires of the same BF    
    rxy = zeros(length(ix),length(ix));
    for jt = 1:length(ix)
        x = lfpDat.LFPseg{ix(jt)}(:,xx==1);%select HITS
        x = x(:);
        for kt = 1:length(ix)
            y = lfpDat.LFPseg{ix(kt)}(:,xx==1);%select HITS
            y = y(:);
            d = corr(x,y);
            rxy(jt,kt) = d;
            if jt==kt
                rxy(jt,kt) = 0;
            end;
        end;
    end;
   [xi,yi] = ind2sub(size(rxy),find(rxy>.7));
   ix = ix(unique([xi yi]));
    
    %% compute median LFP
    x = zeros(length(ix),nsmp,ntrl);
    for jt = 1:length( ix )
        x(jt,:,:) = lfpDat.LFPseg{ix(jt)}(:,xx==1);%select HITS
    end;
    x = squeeze(median(x,1));
    AVG{it} = x;
    
    %% reject trials with dubious activity (>4SDs)
    x2 = AVG{it}(tIx,:);
    M = ones(size(x2,1),1)*mean(x2,1);
    SD = ones(size(x2,1),1)*std(x2,0,1);
    z = (x2-M)./SD;
    z = max(abs(z),[],1);
    
    AVG{it} = AVG{it}(:,z<4);% reject ARTEFACTS
    
    %% apply previous trial selection and rejection to individual microwire data
    ix = find(strcmp(BFlab,BFid(it)));
    for jt = 1:length(ix)
        lfpDat.LFPseg{ix(jt)} = lfpDat.LFPseg{ix(jt)}(:,xx==1);%select HITS
        lfpDat.LFPseg{ix(jt)} = lfpDat.LFPseg{ix(jt)}(:,z<4);% reject ARTEFACTS
    end;
    if size( lfpDat.LFPseg{ix(1)} ) ~= size( AVG{it} )
        error('matrix dimensions must be consistent');
    end;
    fprintf('\n');
end;

%%
figure;
a = zeros(length(AVG),1);
for it = 1:length( AVG )
    [S1,f1] = mtspectrumc( gradient(AVG{it}(tIx,:)')', params3 );
    [S2,f2] = mtspectrumc( gradient(AVG{it}(tIx,:)')', params4 );
    subplot(4,2,it);
    a(it) = gca;
    hold on;
    plot(f1,S1,'b');
    plot(f2,S2,'r');
    lm(it,:) = [min(S1) max(S1) min(S2) max(S2)];
    title(BFid{it});
    xlabel('Frequency (Hz)');
    ylabel('Power (a.u.)');
end;
set(a,'YLim',[min(min(lm)) max(max(lm))]);

%% 
pbins = 0:9:360;
dt = 800;
dt2 = 195;
for it = 33:40%1:length(BFlab) % loop over the BFs       
    
    fprintf([num2str(it),'/',num2str(length(BFlab))]);
       
    [nsmp,ntrl] = size(lfpDat.LFPseg{it});% get dims of input 
    
    %% compute power spectra
    [S1,f1,Serr1] = mtspectrumc( gradient(lfpDat.LFPseg{it}(tIx,:)')', params3);
                    
    [S2,f2,Serr2] = mtspectrumc( gradient(lfpDat.LFPseg{it}(tIx,:)')', params4);
           
    figure;   
    hold on;
    jbfill(f1,(Serr1(1,:)),(Serr1(2,:)),[0 0 .9],[0 0 .9],1,.5);
    jbfill(f2,(Serr2(1,:)),(Serr2(2,:)),[.9 0 0],[.9 0 0],1,.5);
    hold on;
    plot(f1,(S1),'y');
    plot(f2,(S2),'y');
    title(lfpDat.chanLab{it});
    xlabel('Frequency (Hz)');
    ylabel('Power (a.u.)');
    
    %%
    pbins2 = 0:20:360;
    phi = zeros(ntrl*length(tIx),length(pfoi));
    amp = zeros(ntrl*length(tIx),length(afoi));
    for jt = 1:ntrl
        
        fprintf([num2str(jt),'/',num2str(ntrl)]);
        
        x = lfpDat.LFPseg{it}(:,jt)';
        ix = (jt-1)*length(tIx)+1:(jt-1)*length(tIx)+length(tIx);        
        
        d = zeros(length(tIx),length(pfoi));
        parfor kt = 1:length(pfoi)
            ft = filtfilt(fbp{kt},1,[fliplr(x) x fliplr(x)]);
            ft = hilbert(ft);
            ft = ft(nsmp+1:2*nsmp);
            ft = ft(tIx);            
            d(:,kt) = (angle(ft)+pi).*(180/pi);
        end;
        phi(ix,:) = d;
        
        d = zeros(length(tIx),length(afoi));
        for kt = 1:length(afoi)
            ft = filtfilt(fba{kt},1,[fliplr(x) x fliplr(x)]);
            ft = hilbert(ft);
            ft = ft(nsmp+1:2*nsmp);
            ft = ft(tIx);
            d(:,kt) = abs(ft).^2;
        end;
        amp(ix,:) = d;
        
        fprintf('\n');
    end;
        
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
    H = -squeeze(sum(X.*log(X),3));
    n = length(pbins2);
    MI = (log(n)-H)./log(n);                
     
    tic;
    nrand = 1000;
    MI2 = zeros(nrand,length(afoi),length(pfoi));
    parfor nt = 1:nrand
        fprintf([num2str(nt),'/',num2str(nrand)]);
        
        trlp = 1:ntrl;
        ixa = 1:length(tIx);
        % shuffle the amplitude time series against the phase time series
        amp2 = zeros(ntrl*length(tIx),length(afoi));% not this could also be phase instead
        
        f = 0;
        while f<1
            rIx = trlp(randperm(length(trlp)))';% pick random order of trials
            if (sum(diff([trlp' rIx],[],2)==0)==0)
                f = 1;
            end;
        end;
        
        for jt = 1:ntrl        
            
            ix = (rIx(jt)-1)*length(tIx)+1:(rIx(jt)-1)*length(tIx)+length(tIx);
            
            for kt = 1:length(afoi)
                amp2(ixa,kt) = amp(ix,kt);% this could also be phase
            end;
            
            ixa = ixa+length(tIx);
        end;
        
        X = zeros(length(afoi),length(pfoi),length(pbins2));
        for kt = 1:length( pfoi )
            p = phi(:,kt);
            for lt = 1:length( afoi )
                a = amp2(:,lt);
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
    toc;
     
    pval = zeros(size(MI));
    for jt = 1:size(MI,1)
        for kt = 1:size(MI,2)            
            ref = squeeze( MI2(:,jt,kt) );
            pval(jt,kt) = length(find(ref >= MI(jt,kt)))/length(ref);
        end;
    end;
    
    figure;
    pcolor(pfoi,afoi,MI.*(pval<(0.001/length(lfpDat.chanLab))));
    shading interp;lighting phong;
    xlabel('Frequnecy for Phase (Hz)');
    ylabel('Frequency for Amplitude (Hz)');
    cb = colorbar;
    zlab = get(cb,'YLabel');
    set(zlab,'String','MI');
    caxis([0 5e-3]);
    title(lfpDat.chanLab{it});
    
    %%
    % bandpass filter the LFP
    gP = zeros(ntrl,length(tIx));
    gE = zeros(ntrl,length(tIx));
    for jt = 1:ntrl        
        x   = lfpDat.LFPseg{it}(:,jt)';
        x  = filtfilt(b3,1,[fliplr(x) x fliplr(x)]);% filter LFP to obtain gamma-signal       
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
    for jt = 1:ntrl
        
        selIx = findpeaks( -gP(jt,:) );% detect gamma troughs
        selIx = selIx.loc;selIx([1 end]) = [];
        selIx = selIx( gE(jt,selIx)>trsh);%apply treshold
                
        y1 = lfpDat.LFPseg{it}(:,jt)';% raw LFP        
        y2 = filtfilt(b1,1,[fliplr(y1) y1 fliplr(y1)]);% filter LFP to obtain theta-signal
        y3 = filtfilt(b3,1,[fliplr(y1) y1 fliplr(y1)]);% filter LFP to obtain gamma-signal        
        y3 = abs(hilbert(y3)).^2;
        
        y2 = y2(nsmp+1:nsmp*2);        
        y3 = y3(nsmp+1:nsmp*2);
        
        y1 = y1(tIx);
        y2 = y2(tIx);
        y3 = y3(tIx);
         
        % gamma trough triggered LFP
        for kt = 1:length(selIx)
            if (selIx(kt)-dt2 > 0) && (selIx(kt)+dt2<length(x))
                cnt2 = cnt2+1;
                ptA1(cnt2,:) = y1(selIx(kt)-dt2:selIx(kt)+dt2);%raw LFP
                ptA2(cnt2,:) = y2(selIx(kt)-dt2:selIx(kt)+dt2);% theta filtered LFP
                ptA3(cnt2,:) = y3(selIx(kt)-dt2:selIx(kt)+dt2);% gamma filtered LFP
            end;
        end;
        
    end;
    ptA1(cnt2+1:end,:) = [];
    ptA2(cnt2+1:end,:) = [];
    ptA3(cnt2+1:end,:) = [];
    
    x1 = ptA3(find(sign(ptA2(:,find(-dt2:dt2==0)))==-1),find(-dt2:dt2==0));
    x1 = (x1-M)./SD;
    x2 = ptA3(find(sign(ptA2(:,find(-dt2:dt2==0)))==1),find(-dt2:dt2==0));
    x2 = (x2-M)./SD;
    
    ix = find(ptA3(:,find(-dt2:dt2==0)) > trsh);
  
    figure;
    hold on;
    bar([1 2],[mean(x1) mean(x2)]);
    plot([1 1],[mean(x1) mean(x1)+std(x1)/sqrt(length(x1)-1)],'k');
    plot([2 2],[mean(x2) mean(x2)+std(x2)/sqrt(length(x2)-1)],'k');    
    xlim([0 3]);
    set(gca,'XTick',1:2);
    set(gca,'XTickLabel',{'\Theta Trough' '\Theta peak'});
    ylabel('\gamma-power [\sigma]');            
    
    ptA1 = ptA1(ix,:);
    ptA2 = ptA2(ix,:);
    ptA3 = ptA3(ix,:);    
    
    sIx1 = find(sign(ptA3(:,find(-dt2:dt2==0)))==1);
    sIx2 = find(sign(ptA3(:,find(-dt2:dt2==0)))==-1);
    
    figure;
    subplot(4,1,1:3);
    imagesc(-dt2:dt2,1:size(ptA1,1),ptA1([sIx1;sIx2],:));
    axis xy;
    xlabel('Time (ms)');
    ylabel('Frequency (Hz)');
    subplot(4,1,4);
    hold on;
    [ax,h1,h2] = plotyy(-dt2:dt2,mean(ptA1,1),-dt2:dt2,mean(ptA2,1));
    %plot(ax(1),-dt2:dt2,mean(ptA2(sIx1,:),1),'m');
    %plot(ax(1),-dt2:dt2,mean(ptA2(sIx2,:),1),'g');
    axis(ax,'tight');
    axis(ax,'off');    
    
    xc1 = zeros(size(ptA2,1),dt2*2+1);
    xc2 = zeros(size(ptA2,1),dt2*2+1);
    for jt = 1:size(ptA2,1)
        
        [xc1(jt,:),lag] = xcorr(ptA2(jt,:),dt2,'coeff');
        [xc2(jt,:),~] = xcorr(ptA2(jt,:),ptA3(jt,:),dt2,'coeff');
    end;
    M = mean(xc2,1);
    SE = std(xc2,0,1)/sqrt(size(xc2,1)-1);
    
    figure;
    hold on;
    plot(lag,mean(xc1,1),'k');
    jbfill(lag,M-SE,M+SE,'r','r',1,.5);
    hold on;
    plot(lag,mean(xc2,1),'r');
    
    %%
    avg = zeros(ntrl,102,length(pbins));
    %avg = zeros(ntrl,length(afoi),length(pbins));
    ptA4 = zeros(20e3,length(-dt:dt));
    ptA5 = zeros(20e3,length(-dt:dt));
    cnt = 0;
    for jt = 1:ntrl
        fprintf([num2str(jt),'/',num2str(ntrl),'\n']);

        x   = lfpDat.LFPseg{it}(:,jt)';% raw signal
        
        % apply filters
        ft  = filtfilt(b1,1,[fliplr(x) x fliplr(x)]);% theta filtered LFP        
        ft2  = filtfilt(b3,1,[fliplr(x) x fliplr(x)]);% gamma filtered LFP
        
        % hilbert transform -> analytical phase & amp
        ft2 = abs(hilbert(ft2)).^2;% power envelope of gamma filtered LFP
        phi2 = angle(hilbert(ft));% analytical phase of theta filtered LFP
        
        %truncate padding
        ft = ft(nsmp+1:nsmp*2);
        ft2 = ft2(nsmp+1:nsmp*2);        
        phi2 = phi2(nsmp+1:nsmp*2);
        
        % extract samples of interest
        x = x(tIx);
        ft = ft(tIx);
        ft2 = ft2(tIx);
        phi2 = phi2(tIx);
        
        [phii,~,ixT] = interpPhasel(ft);
       
%         phi2 = phi2(ix(1,1):ix(end,end));
%         phi3 = (phi2+pi).*(180/pi);        
        
        % theta trough triggered LFP
        for kt = 1:length(ixT)% loop across troughs
            if (ixT(kt)-dt > 0) && (ixT(kt)+dt < length(x))
                cnt = cnt+1;
                d = x(ixT(kt)-dt:ixT(kt)+dt);
                ptA4(cnt,:) = d;% aligh LFP to theta troughs
                d = ft2(ixT(kt)-dt:ixT(kt)+dt);
                ptA5(cnt,:) = d;% aligh gamma power to theta troughs
            end;
        end;
        
        movingwin2          = [.025 .001];
        T = movingwin2(1);
        W = 40;
        params2.tapers      = [T*W 2*(T*W)-1];
        params2.trialave    = 1;
        params2.pad         = 4;
        [S1,t1,f1] = mtspecgramc( gradient(lfpDat.LFPseg{it}(:,jt)')',movingwin2,params2);
        S1 = S1(t1-5>=-1 & t1-5<=4,:);
        
        %y  = lfpDat.LFPseg{it}(:,jt);
        %for mt = 1:length(afoi)
%             ft = filtfilt(fb{mt},1,[fliplr(y) y fliplr(y)]);
%             ft = abs(hilbert(ft)).^2;
%             ft = ft(nsmp+1:nsmp*2);
%             ft = ft(tIx);
            
            parfor lt = 1:length(pbins)-1
                avg(jt,:,lt) = avg(jt,:,lt)+mean( S1( phii>=pbins(lt) & phii < pbins(lt+1),: ) ,1);
                %avg(jt,mt,lt) = avg(jt,mt,lt)+mean( ft( yi>=pbins(lt) & yi < pbins(lt+1)) );
            end;
            avg(jt,:,end) = avg(jt,:,1);
            %avg(jt,mt,end) = avg(jt,mt,1);
        %end;
        
        fprintf('\n');
    end;
    ptA4(cnt+1:end,:) = [];
    ptA5(cnt+1:end,:) = [];
    
    %%
    M = repmat(mean(mean(avg,3),1),[size(avg,1) 1 size(avg,3)]);
    %x = squeeze(mean(mean(avg(:,f>=50 & f<=100,:),2),3));
    avg = log(avg)-log(M);
    
    x = mean(ptA4,2);
    
    movingwin2          = [.1 .001];
    T = movingwin2(1);
    W = 20;
    params2.tapers      = [T*W 2*(T*W)-1];
    params2.trialave    = 1;
    params2.pad         = 4;
    
    [S2,t2,f2] = mtspecgramc(gradient(ptA4)',movingwin2,params2);
    M = ones(size(S2,1),1)*mean(S2,1);
    S2 = log(S2)-log(M);
    %     [S1,t,f] = mtspecgramc(gradient(ptA4(x<median(x),:))',movingwin2,params2);
    %     M = ones(size(S1,1),1)*mean(S1,1);
    %     S1 = log(S1)-log(M);
    %     [S2,t,f] = mtspecgramc(gradient(ptA4(x>median(x),:))',movingwin2,params2);
    %     M = ones(size(S2,1),1)*mean(S2,1);
    %     S2 = log(S2)-log(M);

    %%
    figure;
    hold on;
    pcolor(pbins,f1,(squeeze(mean(avg,1))));shading interp;   
    %pcolor(pbins,afoi,(squeeze(mean(avg,1))));shading interp;   
    plot([180 180],[min(f1) max(f1)],'w--');axis tight;%caxis([-.5 .5]);
    %plot([180 180],[min(afoi) max(afoi)],'w--');axis tight;%caxis([-.5 .5]);
    set(gca,'XTick',[0:90:360]);                 
    xlabel('\Theta-phase (deg)');
    ylabel('Frequency (Hz)');
    cb = colorbar;
    zlab = get(cb,'YLabel');
    set(zlab,'String','\gamma-power (dB)');
    
    figure;
    subplot(4,1,1:3);
    imagesc(t2.*1e3-dt,f2,S2');axis xy;
    xlim([-300 300]);
    ylabel('Frequency (Hz)');
    subplot(4,1,4);
    [ax,h1,h2] = plotyy(-dt:dt,mean(ptA4,1),t2.*1e3-dt,mean(S2(:,f2>=50 & f2<=100),2));
    set([h1 h2],'LineWidth',3);
    box(ax(1),'off');
    axis(ax,'tight');
    set(ax,'Xlim',[-300 300]);
    caxis([-.1 .1]);
    xlabel('Time (ms)');
    ylabel(ax(1),'LFP-amplitude (\muV)');
    ylabel(ax(2),'\gamma-power (dB)');
    
end;

%%
figure;
x = [];
for it = 1:length(AI)
    x = [x AI{it}'];
end;
m = mean(x);
[n,x] = hist(x,1e2);
figure;
hold on;
bar(x,n,'EdgeColor','b');
plot([m m],[min(n) max(n)],'r--');
xlim([min(x) max(x)]);
ylabel('# of Theta cycles');

%% estimate peaks and troughs based on waveshape
for it = 4%1:length(BFid) % loop over the BFs
    
    [nsmp,ntrl] = size(lfpDat.LFPseg{it});
    phiIx1 = cell(1,ntrl);
    phiIx2 = cell(1,ntrl);
    eP1 = [];   eP2 = [];
    dt = 800;
    for jt = 1:ntrl
        
        x   = lfpDat.LFPseg{it}(:,jt)';
        ft  = filtfilt(b1,1,[fliplr(x) x fliplr(x)]);   
        ft2  = filtfilt(b2,1,[fliplr(x) x fliplr(x)]);
        s   = conv(ft,rectwin(110),'same')./sum(rectwin(110));
        
        ft = ft(nsmp+1:nsmp*2);
        ft2 = ft2(nsmp+1:nsmp*2);
        s = s(nsmp+1:nsmp*2);
        
        x = x(tIx);
        ft = ft(tIx); 
        ft2 = ft2(tIx);
        s = s(tIx);
        
        [ix] = findpeaks( ft );
        [ix2] = findpeaks( ft2 );
        
        sel = zeros(length(ix2.loc),1);
        for kt = 1:length(ix2.loc)            
            [~,mIx] = min(abs(ix2.loc(kt)-ix.loc));
            sel(kt) = mIx;
        end;
        ix.loc = ix.loc(sel);
        ix.loc = ix.loc(sign(ft(ix.loc))==1);
        ix1 = ix.loc;
        
%         ix1 = find(diff(sign(diff(s)))==2)+2;
%         ix2 = find(diff(sign(diff(s)))==-2)+2;
%         ix1 = ix1(sign(s(ix1))==-1);
%         ix2 = ix2(sign(s(ix2))==1);                                       
%         
%         dum = [];
%         dum(ix1) = -1;
%         dum(ix2) = 1;
%         chck = [find(dum~=0)' dum(dum~=0)'];
%         chck = [ chck [2;diff(chck(:,2),[],1)] ];
%         
%         delIx = [];
%         p = find(chck(:,3)==0);
%         while ~isempty(p)
%             ix = min(p);
%             ix = [ix-1 ix];
%             f = 0;
%             while f<1
%                 if (ix(end)+1<=size(chck,1)) && (chck(ix(end)+1,3) == 0)
%                     ix = [ix ix(end)+1];
%                 else
%                     f =1;
%                 end;
%             end;
%             p(find(ismember(p,ix))) = [];
%             
%             [~,mIx] = max(abs(s(chck(ix,1))));
%             del = setdiff(1:length(ix),mIx);
%             delIx = [delIx;ix(del)'];
%         end;
%         chck(delIx,:) = [];
%         ix1 = chck(find(chck(:,2)==-1),1);
%         ix2 = chck(find(chck(:,2)==1),1);
%         
%         if ix2(1)<ix1(1)
%             ix2(1) = [];
%         end;
%         if ix1(end)>ix2(end)
%             ix1(end) = [];
%         end;
        
        %sel = 1:length(ix1);
        %sel = find(diff([ix1 ix2],[],2)>=floor((1/12)/2*1e3) & diff([ix1 ix2],[],2)<=floor((1/3)/2*1e3));
        %ix1 = ix1(sel);
        %ix2 = ix2(sel);
        
%         ix3 =zeros(length(ix1),1);
%         for kt = 1:length(ix1)
%             d = min(find(sign(s(ix1(kt):end))==1));
%             ix3(kt) = ix1(kt)+ d;
%         end;
%         
%         ix4=zeros(length(ix2),1);
%         for kt = 1:length(ix2)
%             d = min(find(sign(s(ix2(kt):end))==-1));
%             ix4(kt) = ix2(kt)+ d;
%         end;
%         
%         x = [ix1(1) ix3(1) ix2(1) ix4(1)];
%         y = [0 90 180 270];
%         xi = ix1(1):ix4(1);
%         yi = interp1(x,y,xi)
%         
%         
%         phiIx1{jt} = ix1;
%         phiIx2{jt} = ix2;
                
        ix = ix1;
        ix = ix(ix-dt>0);
        ix = ix(ix+dt<length(x));
        for lt = 1:length(ix)
            eP1 = [eP1;x(ix(lt)-dt:ix(lt)+dt)];
        end;
        M = mean(eP1,2)*ones(1,size(eP1,2));
        eP1 = eP1-M;
%         
%         ix = ix.loc;
%         ix = ix(ix-dt>0);
%         ix = ix(ix+dt<length(x));
%         for lt = 1:length(ix)
%             eP2 = [eP2;x(ix(lt)-dt:ix(lt)+dt)];
%         end;
%         M = mean(eP2,2)*ones(1,size(eP2,2));
%         eP2 = eP2-M;
%         
%         figure;
%         hold on;
%         plot(ft,'m','LineWidth',3);hold on;
%         plot(s,'k','LineWidth',3);
%         plot(ft2,'g','LineWidth',3);
%         plot(ix1,s(ix1),'bo','MarkerFaceColor','b');
%         plot(ix2,s(ix2),'ro','MarkerFaceColor','r');
%         plot(ix1*ones(1,2),[min(s) max(s)],'k--');
%         plot(ix2*ones(1,2),[min(s) max(s)],'k--');

    end;
        
end;

%%
[Sx,tx,fx] = mtspecgramc( gradient(eP1)', movingwin2, params2 );
%figure;plot(fx,mean(Sx,1));
M = ones(size(Sx,1),1)*mean(Sx,1);
Sx = log(Sx)-log(M);

figure;
subplot(4,1,[1:3]);
imagesc(tx-dt./1e3,fx,Sx');
axis xy;
set(gca,'Xlim',[-0.5 0.5]);
subplot(4,1,4);
hold on;
[ax,h1,h2] = plotyy(-dt:dt,mean(eP1,1),tx-dt./1e3,mean(Sx(:,fx>=60 & fx<=80),2));
axis(ax,'tight');
set(ax(2),'Xlim',[-0.5 0.5]);
set(ax(1),'Xlim',[-500 500]);

%%
[Sx,tx,fx] = mtspecgramc( gradient(eP2)', movingwin2, params2 );
%figure;plot(fx,mean(Sx,1));
M = ones(size(Sx,1),1)*mean(Sx,1);
Sx = log(Sx)-log(M);

figure;
subplot(4,1,[1:3]);
imagesc(tx-dt./1e3,fx,Sx');
axis xy;
set(gca,'Xlim',[-0.5 0.5]);
subplot(4,1,4);
hold on;
[ax,h1,h2] = plotyy(-dt:dt,mean(eP2,1),tx-dt./1e3,mean(Sx(:,fx>=60 & fx<=80),2));
set(ax(2),'Xlim',[-0.5 0.5]);
set(ax(1),'Xlim',[-500 500]);


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
for zt = 26%1:length( LFPseg ) 
           
    [AVG] = lfpDat.LFPseg{zt};
    
    %%
    ntrl = size(AVG,2);
    
    missIdx2 = find(ismember(sort([lfpDat.missIdx;lfpDat.hitIdx]),lfpDat.missIdx));
    hitIdx2 = find(ismember(sort([lfpDat.missIdx;lfpDat.hitIdx]),lfpDat.hitIdx));
    xx = zeros(ntrl,1);
    xx(hitIdx2) = 1;
    
    %%    
    AVG = AVG(:,xx == 1);
    M = ones(size(AVG,1),1) * mean(AVG,1);
    SD = ones(size(AVG,1),1) * std(AVG,0,1);
    AVG = (AVG-M)./SD;   
    
    z = max(abs(AVG),[],1);    
    AVG = AVG(:,z<4);
    
    %%
    [trl] = spkDat.sortedSpikes{zt}.trl;  
    [ts]  = spkDat.sortedSpikes{zt}.SpikeTimesSeg.*1e3;
    
    selIx = find( ismember(trl,find(xx==1)) );
    trl = trl(selIx);
    ts  = ts(selIx);
    trl = trl(ts>= -1e3 & ts <= 4e3 );
    ts = ts(ts>= -1e3 & ts <= 4e3 );
    
    dt = -1e3:250:4e3;
    spkC = zeros(ntrl,length(dt));
    for kt = 1:ntrl
        ix = find(trl == kt);
        x = ts(ix);
        spkC(kt,:) = hist(x,dt);
    end;
    spkC = spkC(z<4,:);
    
    
%     %%
%     x = AVG;
%     M = ones(size(x',1),1)*mean(x',1);
%     x = x'-M;
%     
%     x = gradient( x )';
%     
%     [s1,t1,f1] = mtspecgramc( x, movingwin1 , params1);
%     [s2,t2,f2] = mtspecgramc( x , movingwin2, params2  );
%     t1 = t1-5; t2 = t2-5;
%     s1 = s1(t1>=-1 & t1 <=4,:);
%     s2 = s2(t2>=-1 & t2 <=4,:);
%     t1 = t1(t1>=-1 & t1 <=4);
%     t2 = t2(t2>=-1& t2 <=4);
%     
%     %%
%     x = AVG(trlTime>=-1 & trlTime <=4,:);
%     M = mean(x,2)*ones(1,size(x,2));
%     
%     [s3,f3] = mtspectrumc( x, params3  );
%     [s4,f4] = mtspectrumc( x , params4 );
%     s3 = s3(f3>=3);
%     s4 = s4(f4>=3);
%     f3 = f3(f3>=3);
%     f4 = f4(f4>=3);
%     
%     b = regress(log10(s3),[ones(size(f3')) log10(f3')]);
%     yp = [ones(size(f3')) log10(f3')]*b;
%     s3c = 10.^(log10(s3)-yp);
%     
%     b = regress(log10(s4),[ones(size(f4')) log10(f4')]);
%     yp = [ones(size(f4')) log10(f4')]*b;
%     s4c = 10.^(log10(s4)-yp);
%     
%     %%
%     figure;
%     subplot(4,1,1:3);
%     imagesc(t2,f2,s2');axis xy;
%     title(['n= ',num2str(size(AVG,2))]);
%     subplot(4,1,4);
%     imagesc(t1,f1,s1');axis xy;
%         
%     figure;
%     subplot(211);
%     hold on;
%     plot(f3,20*log10(s3));
%     plot(f4,20*log10(s4),'r');
%     xlim([0 100]);
%     subplot(212);
%     hold on;
%     plot(f3,s3c);
%     plot(f4,s4c,'r');
%     xlim([0 100]);
    
    %%
    bpf = [4 11];
    [b1] = fir1(3*floor(lfpDat.Fs/bpf(1)),[ bpf ]./(lfpDat.Fs/2),'bandpass');
    
    bpf = [1 30];
    [b2] = fir1(3*floor(lfpDat.Fs/bpf(1)),[ bpf ]./(lfpDat.Fs/2),'bandpass');
    
    dt = 750;
    c = 0;
    c2 = 0;
    eP = [];
    eP2 = [];
    for kt = 1:size(AVG,2);
                 
        x = AVG(:,kt);
        x = x(lfpDat.trlTime >=-1 & lfpDat.trlTime <=4);        
        n = length(AVG(:,kt));
        
        d = filtfilt(b2,1,[fliplr(AVG(:,kt))' AVG(:,kt)' fliplr(AVG(:,kt))']);
        d = d(n+1:2*n);        
        d = d(lfpDat.trlTime >=-1 & lfpDat.trlTime <=4);                       
                
        [ix] = findpeaks(d);
        [ix3] = findpeaks(-d);
        
        n = length(AVG(:,kt));
        d2 = filtfilt(b1,1,[fliplr(AVG(:,kt))' AVG(:,kt)' fliplr(AVG(:,kt))']);
        d2 = d2(n+1:2*n);
        d2 = d2(lfpDat.trlTime >=-1 & lfpDat.trlTime <=4);
        [ix2] = findpeaks(d2);
        [ix4] = findpeaks(-d2);
        ix2.loc([1 end]) = [];
        ix4.loc([1 end]) = [];
        
        sel = [];
        for lt = 1:length(ix2.loc)
            [~,mIx] = min(abs(ix2.loc(lt)-ix.loc));
            sel(lt) = mIx;
        end;
        ix.loc = ix.loc(sel);
        
        sel = [];
        for lt = 1:length(ix4.loc)
            [~,mIx] = min(abs(ix4.loc(lt)-ix3.loc));
            sel(lt) = mIx;
        end;
        ix3.loc = ix3.loc(sel);                
        
        while ix.loc(1)<ix3.loc(1)
            ix.loc(1) = [];
        end;
        
         while ix3.loc(end)>ix.loc(end)
            ix3.loc(end) = [];
        end;
                
        
%         figure;plot(d);hold on;plot(d2,'r');
%         plot(ix3.loc,d2(ix3.loc),'go');
%         plot(ix4,d2(ix4),'ko');
%         plot(ix.loc,d2(ix.loc),'mo');
%         plot(ix2,d2(ix2),'co');
        
                                
        for lt = 1:length(ix.loc)
            if (ix.loc(lt)-dt>0) && (ix.loc(lt)+dt<length(x))
                c = c+1;
                eP(c,:) = d(ix.loc(lt)-dt:ix.loc(lt)+dt);
            end;
        end;
        
        for lt = 1:length(ix3.loc)
            if (ix3.loc(lt)-dt>0) && (ix3.loc(lt)+dt<length(x))
                c2 = c2+1;
                eP2(c,:) = d(ix3.loc(lt)-dt:ix3.loc(lt)+dt);
            end;
        end;
        
    end;
    
    %%
    bpf = [4 8];
    [b2] = fir1(3*floor(lfpDat.Fs/bpf(1)),[ bpf ]./(lfpDat.Fs/2),'bandpass');
    
    c2 = 0;
    c3 = 0;
    eP = [];
    FR = [];
    for it = 1:size(AVG,2)
        
        ft = AVG(lfpDat.trlTime >=-1 & lfpDat.trlTime <=4,it);

        n = length(AVG(:,it));
        d = filtfilt(b2,1,[fliplr(AVG(:,it))' AVG(:,it)' fliplr(AVG(:,it))']);
        %d = angle(hilbert(d));
        d = d(n+1:2*n);
                
        ft2 = d(lfpDat.trlTime >=-1 & lfpDat.trlTime <=4);
        [ix] = findpeaks(ft2);
        %ix.loc = ix.loc(sign(ft2(ix.loc))==1);
        %ix.loc = ix.loc(ft2(ix.loc) < -pi+0.01);
        
%         c = 0;
%         ix2 = 1:round(1/bpf(2).*1e3);
%         n = floor(length(ft2)/round(1/bpf(2).*1e3));
%         sel = [];
%         for jt = 1:n            
%             ix3 = find(ismember(ix.loc,ix2));  
%             if c>1
%                 ix3 = ix3(abs(sel(c)-ix.loc(ix3))>round(1/bpf(2).*1e3));
%             end;
%             if ~isempty(ix3)
%                 c = c+1;
%                 [~,d] = max(ft2(ix.loc(ix3)));
%                 sel(c) = ix.loc(ix3(d));
%             end;
%             ix2 = ix2+round(1/bpf(2).*1e3);            
%         end;
        sel = ix.loc;
        
        dt = 800;
        for kt = 1:length(sel)
            if (sel(kt) -dt>0) && (sel(kt)+dt<length(ft))
                c2 = c2+1;
                eP(c2,:) = ft(sel(kt)-dt:sel(kt)+dt);
            end;
        end;        
        
    end;
    M = mean(eP,2)*ones(1,size(eP,2));
    eP = eP-M;
%     p = mean(abs( eP ).^2,2);
%     lm = prctile(p,95);
%     eP = eP(p>lm,:);

    [Sx,tx,fx] = mtspecgramc( gradient(eP1)', movingwin2, params2 );
    %figure;plot(fx,mean(Sx,1));
    M = ones(size(Sx,1),1)*mean(Sx,1);
    Sx = log(Sx)-log(M);
    
    figure;
    subplot(4,1,[1:3]);
    imagesc(tx-dt./1e3,fx,Sx');
    axis xy;
    set(gca,'Xlim',[-0.5 0.5]);
    subplot(4,1,4);
    hold on;
    [ax,h1,h2] = plotyy(-dt:dt,mean(eP1,1),tx-dt./1e3,mean(Sx(:,fx>=60 & fx<=80),2));
    set(ax(2),'Xlim',[-0.5 0.5]);
    set(ax(1),'Xlim',[-500 500]);
    
    %%
    bpf = [60 80];
    [b2] = fir1(3*floor(lfpDat.Fs/bpf(1)),bpf./(lfpDat.Fs/2),'bandpass');
    
    c2 = 0;
    c3 = 0;
    eP = [];
    FR = [];
    for it = 1:size(AVG,2)
        
        n = length(AVG(:,it));
        
        ft = AVG(:,it);
        ft = ft(lfpDat.trlTime >=-1 & lfpDat.trlTime <=4);
        
        d = filtfilt(b2,1,[fliplr(AVG(:,it))' AVG(:,it)' fliplr(AVG(:,it))']);
        %d = angle(hilbert(d));
        %d = d.*(180/pi);
        %d = d+180;
        d = abs(hilbert(d)).^2;
        d = d(n+1:2*n);
                
        ft2 = d(lfpDat.trlTime >=-1 & lfpDat.trlTime <=4);
        
        m = mean(abs(ft2));
        sd = 3*std(abs(ft2));
        
        trsh = abs(ft2) > m+sd;
        seg1 = find(diff(trsh)==1)+1;
        seg2 = find(diff(trsh)==-1)+1;
        delIx = [];
        c = 0;
        for kt = 1:length(seg2)
            if sum(seg2(kt)<seg1) == length( seg1 )
                c = c+1;
                delIx(c) = kt;
            end;
        end;
        seg2(delIx) = [];
        
        ix = [];
        ix.loc =[];
        for kt = 1:size(seg,1)
            rg = seg(kt,1):seg(kt,2);
            [~,mx] = max( ft2(rg) );
            ix.loc = [ix.loc rg(mx)];
        end;
        
        %[ix.loc] = [find(d >=265) find(d <=275)];%
        %[ix] = findpeaks(ft2.^2);
        
%         c = 0;
%         ix2 = 1:125;
%         n = floor(length(ft2)/125);
%         sel = [];
%         for jt = 1:n            
%             ix3 = find(ismember(ix.loc,ix2));  
%             if c>1
%                 ix3 = ix3(abs(sel(c)-ix.loc(ix3))>125);
%             end;
%             if ~isempty(ix3)
%                 c = c+1;
%                 [~,d] = max(ft2(ix.loc(ix3)));
%                 sel(c) = ix.loc(ix3(d));
%             end;
%             ix2 = ix2+125;            
%         end;
        sel = ix.loc;
        dt = 400;
        for kt = 1:length(sel)
            if (sel(kt) -dt>0) && (sel(kt)+dt<length(ft))
                c2 = c2+1;
                eP(c2,:) = ft(sel(kt)-dt:sel(kt)+dt);
            end;
        end;
        
        for kt = 1:length(sel)
            if (sel(kt) -dt>0) && (sel(kt)+dt<size(spkC,2))
                c3 = c3+1;
                FR(c3,:) = spkC(it,sel(kt)-dt:sel(kt)+dt);
            end;
        end;
        
    end;
    p = mean(eP.^2,2);
    lm = prctile(p,[25 75]);
    sel = find(p >= lm(1) & p <= lm(2));
    
    bpf = [15];
    [b2] = fir1(3*floor(lfpDat.Fs/bpf(1)),bpf./(lfpDat.Fs/2),'low');
    x = mean(eP(sel,:),1);
    n = length(x);
    d = filtfilt(b2,1,[fliplr(x) x fliplr(x)]);
    d = d(n+1:2*n);
    eP2= d;
    
    [~, sIx] = sort(eP(sel,dt+1));
    figure;
    subplot(221);
    imagesc(-dt:dt,1:size(eP(sel,:),1),eP(sel(sIx),:));
    axis xy;
    subplot(223);
    hold on;
    plot(-dt:dt,mean(eP(sel,:),1));
    plot(-dt:dt,mean(eP2,1),'r');
    xlim([-dt dt]);
    
    x = (eP(sel,:))';
    [Sx1,fx1] = mtspectrumc(x,params1);
    [Sx2,fx2] = mtspectrumc(x,params2);
    subplot(222);
    plot(fx1,20*log10(mean(Sx1,2)));hold on;plot(fx2,20*log10(mean(Sx2,2)),'r');
    xlim([0 100]);
    
    x = gradient(eP)';
    [Sx1,fx1] = mtspectrumc(x,params1);
    [Sx2,fx2] = mtspectrumc(x,params2);
    subplot(224);
    plot(fx1,mean(Sx1,2));hold on;plot(fx2,mean(Sx2,2),'r');
    xlim([0 100]);
    
end;

%%
x = AVG;
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