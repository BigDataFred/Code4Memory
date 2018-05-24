function [thetaTmstmps] = detectThetaBouts(LFPsig,Fs,thetaFreq)
%%
%load('P09_fVSpEM_2017-08-28_14-20-18_medianANDicaLFP.mat','AVG');
%load('P09_fVSpEM_2017-08-28_14-20-18_medianLFP_powerSpectra.mat');

%%
Fs = Fs;
bpf = [thetaFreq(1) thetaFreq(2)];
[b] = fir1(3*(Fs/bpf(1)),bpf./(Fs/2),'bandpass');

%%
[trlTime] = -5:1/Fs:7;
tIx = find(trlTime>=-1 & trlTime <=4);

%%
movingwin = [.8 .2];

T = movingwin(1);
W = 2;
k = floor(2*(T*W)-1);

params                  = [];
params.pad              = 2;
params.Fs               = Fs;
params.fpass            = [0 30];
params.tapers           = [T*W k];
params.trialave         = 1;

%%
[nsmp,~] = size( LFPsig );
[thetaTmstmps] = cell( size(LFPsig,2) ,1);

%%
%figure;
for jt = 1:size(LFPsig,2)
    
    x = LFPsig(:,jt);
    
    x=locdetrend(x,Fs,[1 .25]);
    x = x-mean(x);
    
    %[AR,RC,PE] = mvar(x, 2, 1);
    
    ft = filtfilt(b,1,[fliplr(x) x fliplr(x)]);
    ft = ft(nsmp+1:2*nsmp);
    
    [S,t,f] = mtspecgramc( gradient(x')', movingwin, params);
    t = t-5;
    S = S(t>=-1&t<=4,:);
    t = t(t>=-1&t<=4);
    pow1 = mean(S(:,f>=bpf(1) & f<=bpf(2)),2);
    pow2 = mean(S(:,f>=1 & f<=3),2);
    pow3 =mean(S(:,f>=12&f<=20),2);
    
    pr1 = pow1./pow2;
    pr2 = pow1./pow3;
    d1 = pr1>2;
    d2 = pr2>2;
    trsh = (d1+d2) >=1;
    dx = [0 diff(trsh)'];
    ix1 = find(sign(dx)==1);
    ix2 = find(sign(dx)==-1);
    
    ix = [];
    if ~isempty(ix1) && ~isempty(ix2)
        if length(ix2)>length(ix1)
            ix1 = [1 ix1];
        end;
        if (min(ix1)>min(ix2)) && (length(ix1)==length(ix2))
            ix1 = [1 ix1];
            ix2 = [ix2 length(t)];
        end;
        if length(ix1)>length(ix2)
            ix2 = [ix2 length(t)];
        end;
        
        ix = [ix1' ix2'];
        ix = ix(diff(t(ix),[],2)/(1/mean(bpf))>2,:);
    end;
    
%     subplot(6,1,1:3);
%     a = gca;
%     hold on;
%     imagesc(t,f,(S)');axis xy;
%     for lt = 1:size(ix,1)
%         plot(t(ix(lt,1))-.1*ones(1,2),[min(f) max(f)],'b');
%         plot(t(ix(lt,2))+.1*ones(1,2),[min(f) max(f)],'r');
%     end;
%     axis tight;
%     set(gca,'XLim',[min(t) max(t)]);
%     subplot(6,1,4);
%     hold on;
%     plot(trlTime(tIx),(x(tIx)));
%     plot(trlTime(tIx),ft(tIx),'g','LineWidth',3);
%     for kt = 1:size(ix,1)
%         h = area([t(ix(kt,1))-.1 t(ix(kt,2))+.1],max(x(tIx))*ones(1,2),min(x(tIx)));
%         set(h,'FaceColor',[.8 .8 .8],'FaceAlpha',.4);
%     end;
%     axis(gca,'tight');axis(gca,'tight');set(gca,'XLim',[min(t) max(t)]);
%     set(gca,'Xtick',get(gca,'XTick'));
%     subplot(6,1,5);
%     [ax,h1,h2] = plotyy(t,dx,t,mean(S(:,f>=bpf(1) & f<=bpf(2)),2).*trsh);
%     hold(ax(1),'on');
%     set(ax(1),'Ylim',[-1.5 1.5]);
%     if ~isempty(ix)
%         plot(ax(1),t(ix(:,1)),dx(ix(:,1)),'bo');
%         plot(ax(1),t(ix(:,2)),dx(ix(:,2)),'ro');
%     end;
%     set(ax,'Xtick',get(a,'XTick'));
%     set(ax,'XLim',[min(t) max(t)]);
%     subplot(6,1,6);
%     hold on;
%     plot(t,pow1./pow2,'b');
%     plot(t,pow1./pow3,'r');
%     set(gca,'Xtick',get(a,'XTick'));
%     set(gca,'XLim',[min(t) max(t)]);
    
    tTs = [];
    if ~isempty(ix)
        tTs =  t(ix);
        tTs(:,1) = tTs(:,1)-.1;
        tTs(:,2) = tTs(:,2)+.1;
        ix2 = 1:size(tTs,1)-1;
        cnt = 0;
        while ~isempty(ix2)
            cnt = cnt+1;
            %forward merging
            if round(tTs(cnt,2)*1e3)/1e3>=round(tTs(cnt+1,1)*1e3)/1e3
                x1 = [tTs(cnt,1) tTs(cnt+1,2)];
                tTs(cnt+1,:) = [];             
                tTs(cnt,:) = x1;
            end;
            %backward merging
            if (cnt>1) && (round(tTs(cnt-1,2)*1e3)/1e3>=round(tTs(cnt,1)*1e3)/1e3)
                x1 = [tTs(cnt-1,1) tTs(cnt,2)];
                tTs(cnt,:) = [];
                tTs(cnt-1,:) = x1;
            end;
            ix2(ismember(ix2,[cnt cnt+1])) = [];
        end;
    end;
    
    thetaTmstmps{jt} = tTs;
%     pause;
%     clf;
end;
