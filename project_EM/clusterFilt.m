function [spkDat,BFinf,glFR,FR,uSelIx,instFR,wvfStats,spkRaster] = clusterFilt(x,spkDat,uIx,trlENC,seshLab,plt)
%function [spkDat,BFinf,glFR,FR,uSelIx,trXC,instFR,wvfStats] = clusterFilt(x,spkDat,uIx,trlENC,seshLab,plt)

t = [-2 5]; % time window of interest
t = t.*1e3; % convert s to ms
dt = t(1):250:t(2);% bins used for the isi range
ct = 0;% init counter
uSelIx = [];% init selection index
glFR = {};
FR = {};
trXC = {};
instFR = {};

%%
chLab = cell(1,size(x,1));
wvfStats = cell(1,size(x,1));
for it = 1:size( x ,1 ) % loop over n = MW*units
    fprintf([num2str(it),'/',num2str(size( x ,1 ))]);
    [selDat] = extractClusterDat(spkDat,uIx(it,:));% extracts data of interest
    
    if unique(selDat.unit{1}) ~=0
        ts = selDat.timestamp{1};% raw spiketimes
        ts3 = selDat.time{1}.*1e3;
        trl = selDat.trial{1}; % trial labels
        
        % save the MW label
        chLab{it} =selDat.hdr.label;
        chLab{it}(regexp(chLab{it},'_')) = [];
        
        if ~isempty(ts)
            ts =ts./1e3; % convert raw ts from us to ms
            
            % spike time autocorrelation - returns coincidences
            lag = [];   xc = zeros(length(trlENC),501);            
            ts2 = ts; % retain spike times corresponding to selected trial
            if ~isempty(ts2) && length(ts2)>1
                ts2 = ts2-ts2(1); % normalize time range
                tt =min(ts2):max(ts2);% create the bins
                xn = hist(ts2,tt);% spike count over bins
                [dum] = xcorr(xn,250);% compute autocorr of spike times
                xc = dum;
            end;
            lag = linspace(-250,250,size(xc,2));
%FIXME            
%             nfft = 2^nextpow2(length(xc));
%             y = fft(conv(xc,hanning(nfft),'same'))./1e3;
%             y = y.*conj(y);
%             y = fftshift(y);
%             y = y(nfft/2:end);
%             f = 500*linspace(0,1,nfft/2+1);
%             [~,xcIx] = max(xc);
%             if f(xcIx) >=49 && f(xcIx)<51
%                 return;
%             end;
            
            % adjust autocorr for 0-lag
            xc(max(lag)+1) = NaN;
            % only keep posivite lags
            ix = (length(lag)-1)/2+1:length(lag);
            xc = xc(ix);
            lag = lag(ix);                      
            
            % compute the ISI
            isi = diff(ts);
            pct = length(find(isi<=3))/length(isi);
            isi(isi>250) = [];
            [nIsi,~] = hist(isi,0:250);
            % compute the % of isis below 3 ms
            
            pct = pct.*100;
            
            % compute the instataneous firing rate
            n = [];
            parfor jt = 1:size( x ,2 )
                dt;
                ts = x{it,jt};
                ts(ts<dt(1)) = [];
                ts(ts>dt(end)) = [];
                [n(jt,:),~] = hist( ts  , dt );
            end;
            
            tInt = [-0.95 -0.05;...
                0.3   1.2;...
                2.5   3.5].*1e3;
            
            fr = zeros(length(trlENC),size(tInt,1));
            glfr = zeros(length(trlENC),1);
            parfor jt = 1:length(trlENC)
                ts3;
                dum = ts3(trl == trlENC(jt));
                glfr(jt) = length(find(dum>=min(min(tInt)) & dum <= max(max(tInt))))/(diff([min(min(tInt)) max(max(tInt))])./1e3);
                dum2 = zeros(1,size(tInt,1));
                for kt = 1:size(tInt,1)
                    T = abs(diff(tInt(kt,:)))./1e3;
                    dum2(kt) = length(find(dum>=tInt(kt,1) & dum < tInt(kt,2)))/T;
                end;
                fr(jt,:) = dum2;
            end;
            
            [h1,pval1] = ttest(fr(:,1),fr(:,2));
            [h2,pval2] = ttest(fr(:,1),fr(:,3));
            h1(isnan(h1)) = 0;
            h2(isnan(h2)) = 0;
            
            [snr] = max(sqrt(mean(squeeze(selDat.waveform{1}),2).^2))./mean(selDat.std(:,1));
            [wvfStats{it}] = computeWVFstats( squeeze(selDat.waveform{1}) );
            wvfStats{it}.snr2 = snr;
            wvfStats{it}.nIsi = nIsi;
            wvfStats{it}.pct = pct;
            wvfStats{it}.lag = lag;
            wvfStats{it}.xc = xc;
            
            y = sum(n,1)/size(n,1)/250e-3;
            [~,mIx] = max(y);
                        
            % aply some selection criteria to filter spurious from actual units
            % (pct<=3) && 
            
            if (snr>2 && wvfStats{it}.snr > 0.7) && (dt(mIx) ~=0) && (wvfStats{it}.peakDur > 0.35 && ~isnan(wvfStats{it}.peakDur)) && ( wvfStats{it}.rpct > 0.5)
                
                ct = ct+1;
                uSelIx(ct,1) = it; % memorize non-spurious units
                uSelIx(ct,2) = pval1;
                uSelIx(ct,3) = pval2;
                
                glFR{ct} = glfr;
                FR{ct} = fr;
                instFR{ct} = y;
                
%                 if median(glfr) >8
%                     tt = [selDat.trialtime(1,1) selDat.trialtime(1,2)].*1e3;
%                     tt = tt(1):tt(2);
%                     [trXC{ct}.xc,trXC{ct}.fr,trXC{ct}.lag,trXC{ct}.ntrl] = time_resolved_xcorr(tt,480,trlENC,selDat,[]);
%                 else
%                     trXC{ct}.xc = {};
%                     trXC{ct}.fr = {};
%                     trXC{ct}.lag = {};
%                     trXC{ct}.ntrl = {};
%                 end;
                
                if strcmp(plt,'y') || strcmp(plt,'yes')
                    figure;
                    subplot(4,5,[4 5 9 10]);
                    hold on;
                    plot(spkDat{1}.waveformtime,squeeze(selDat.waveform{1}),'r');
                    plot(spkDat{1}.waveformtime,squeeze(mean(selDat.waveform{1},3)),'k');
                    axis tight;
                    title([chLab{it},',C#',num2str(uIx(it,2))]);
                    ylabel('Amplitude [\muV]');
                    xlabel('Time [ms]');
                    
                    subplot(4,5,[1 2 6 7 11 12]);
                    hold on;
                    for jt = 1:size( x ,2 )
                        r = x{it,jt};
                        t = jt*ones(1,length(r));
                        r = [r;r];
                        t = [t-.5;t+.5];
                        line(r,t,'Color','k');
                        r = [];
                    end;
                    plot([0 0],[0 size(x,2)+1],'Color','r');
                    plot([2e3 2e3],[0 size(x,2)+1],'Color','r');
                    xlim([-.5 5].*1e3);
                    ylim([0 size(x,2)+1]);
                    title([seshLab]);
                    ylabel('Trial #');
                    xlabel('Time [ms]');
                    
                    subplot(4,5,[14 15]);
                    bar(0:250,nIsi);
                    axis tight;
                    title(['nSpikes: ',num2str(sum(sum(n))),', pct(<3ms): ',num2str(round(pct*100)/100)]);
                    ylabel('Count');
                    xlabel('ISI [ms]');
                    
                    subplot(4,5,[19 20]);
                    plot(lag,xc);
                    axis tight;
                    ylabel('Coincidences');
                    xlabel('Lag [ms]');
                    
                    subplot(4,5,[16 17]);
                    y = sum(n,1)/size(n,1)/250e-3;
                    hold on;
                    plot(dt(1:end-1),y(1:end-1),'k-s','LineWidth',3);
                    plot([0 0],[0 max(y)+1],'Color','r');
                    plot([2e3 2e3],[0 max(y)+1],'Color','r');
                    axis tight;xlim([-.5 5].*1e3);
                    ylabel('Spikes/s');
                    xlabel('Time [ms]');
                    
%                     if median(trXC{ct}.fr) >= 5
%                         figure;
%                         hold on;
%                         imagesc(tt,trXC{ct}.lag,trXC{ct}.xc);
%                         plot([0 0],[min(trXC{ct}.lag) max(trXC{ct}.lag)],'w');
%                         plot([2e3 2e3],[min(trXC{ct}.lag) max(trXC{ct}.lag)],'w');
%                         axis tight;
%                         xlabel('Time [ms]');
%                         ylabel('Lag [ms]');
%                         title(['Median firing rate:',num2str(median(trXC{ct}.fr)),'Hz']);
%                     end;
                    
                end;
            end;
        end;
    end;
    fprintf('\n');
end;

%%
[empT1,empT2,pval1,pval2,thr1,thr2] = extractSigRespsFR(FR,1e4);% empT1 = cue vs base, empT2 = encoding vs base
uSelIx(:,4:9) = [empT1' empT2' pval1 pval2 thr1 thr2];

%% only keep data for putative non-spurious units
chLab = chLab(uSelIx(:,1)); % MW label of each unit
wvfStats = wvfStats(uSelIx(:,1));
spkRaster = x(uSelIx(:,1),:);

uIx = uIx(uSelIx(:,1),:);% index of each unit

uID = unique(uIx(:,1));
ix = cell(1,length(uID));
parfor it = 1:length(uID)
    uIx;
    ix{it} = uIx(uIx(:,1) == uID(it),2);
end;

%%
% only keep the MW and clusters that are flagged as non-spurious
BFinf = {};
spkDat2 = {};
chLab2 = {};
cnt = 0;
for it = 1:length(uID)
    
    for jt = 1:length(ix{it})
        cnt = cnt+1;
        spkDat2{cnt}.label                = spkDat{uID(it)}.label([ix{it}(jt)]);
        spkDat2{cnt}.timestamp            = spkDat{uID(it)}.timestamp([ix{it}(jt)]);
        spkDat2{cnt}.waveform             = spkDat{uID(it)}.waveform([ix{it}(jt)]);
        spkDat2{cnt}.unit                 = spkDat{uID(it)}.unit([ix{it}(jt)]);
        spkDat2{cnt}.time                 = spkDat{uID(it)}.time([ix{it}(jt)]);
        spkDat2{cnt}.trial                = spkDat{uID(it)}.trial([ix{it}(jt)]);
        spkDat2{cnt}.hdr                  = spkDat{uID(it)}.hdr;
        spkDat2{cnt}.params               = spkDat{uID(it)}.params;
        spkDat2{cnt}.dimord               =spkDat{uID(it)}.dimord;
        spkDat2{cnt}.cfg                  =spkDat{uID(it)}.cfg;
        spkDat2{cnt}.waveformtime         =spkDat{uID(it)}.waveformtime;
        spkDat2{cnt}.std                  =spkDat{uID(it)}.std;
        spkDat2{cnt}.trialtime            =spkDat{uID(it)}.trialtime;
        
        BFinf(1,cnt) = {spkDat{uID(it)}.hdr.label};
        x = str2double(BFinf{1,cnt}(regexp(BFinf{1,cnt},'\d{1}')));
        BFinf(2,cnt) = {x};
        
        chLab2(cnt) = {spkDat{uID(it)}.hdr.label};
    end;    
end;
chLab = chLab2;

%%
% for it = 1:length(uID)
%     
%     spkDat{uID(it)}.label                = spkDat{uID(it)}.label([ix{it}]);
%     spkDat{uID(it)}.timestamp            = spkDat{uID(it)}.timestamp([ix{it}]);
%     spkDat{uID(it)}.waveform             = spkDat{uID(it)}.waveform([ix{it}]);
%     spkDat{uID(it)}.unit                 = spkDat{uID(it)}.unit([ix{it}]);
%     spkDat{uID(it)}.time                 = spkDat{uID(it)}.time([ix{it}]);
%     spkDat{uID(it)}.trial                = spkDat{uID(it)}.trial([ix{it}]);
%     
%     BFinf(1,it) = {spkDat{uID(it)}.hdr.label};
%     x = str2double(BFinf{1,it}(regexp(BFinf{1,it},'\d{1}')));
%     BFinf(2,it) = {x};
%     
% end;

%%
%[spkDat] = spkDat(uID);
[spkDat] = spkDat2;
if length(chLab) ~= length(spkDat)
    error('events must correspond');
end;
if length(BFinf(1,:)) ~= length(spkDat)
    error('events must correspond');
end;

%%
if strcmp(plt,'y') || strcmp(plt,'yes')
    figure;
    n = length(FR);
    k = round(length(FR)/2);
    ct = 0;
    for it = 1:length(FR)
        ct =ct+1;
        subplot(k,2,it);
        hold on;
        bar([1 2 3 ],mean(FR{it},1),'w');
        for jt = 1:size(FR{it},2)
            M = mean(FR{it}(:,jt));
            SD = std(FR{it}(:,jt))./sqrt(size(FR{it},1)-1);
            plot([jt jt],[M-SD M+SD],'k','LineWidth',2);
            plot([jt-.1 jt+.1],[M+SD M+SD],'k','LineWidth',2);
            plot([jt-.1 jt+.1],[M-SD M-SD],'k','LineWidth',2);
        end;
        ylabel('Firing rate (Hz)');
        set(gca,'XTick',1:3);
        set(gca,'XTickLabel',{'Bas','Cue','Enc'});
        set(gca,'XTickLabelRotation',-35);
        title(chLab{it});
    end;
end;
