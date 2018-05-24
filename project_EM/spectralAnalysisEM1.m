function spectralAnalysisEM1(pID,expMode)

%%
if nargin == 0
    pID = 'P09';%dataset label
    expMode = 'fVSpEM';
end;

%%
addpath(genpath('/home/rouxf/tbx/chronux_2_11/'));
addpath(genpath('/home/rouxf/tbx/eeglab14_1_1b/functions/sigprocfunc/'));
addpath(genpath('/home/rouxf/prj/Bham/code/mcode/utils/'));

%% parameters for spectral analysis (I)
[Fs] = 1e3;

%lower frequencies
%movingwin1 = [1 0.01];
movingwin1 = [.8 .2];
T = movingwin1(1);% length of time window in s
W = 2;% smoothing (+/- 1Hz)
TW = T*W;% time-bandwidth product
k = floor(2*(T*W)-1);% number of tapers

params1              = [];
params1.Fs           = Fs;
params1.pad          = 2;
params1.fpass        = [0.5 30];
params1.tapers       = [TW k];
params1.trialave     = 1;

%higher frequencies
%movingwin2 = [0.25 0.01];
movingwin2 = [0.25 0.0625];
T = movingwin2(1);% length of time window in s
W = 8;% smoothing (+/- 10Hz)cd
TW = T*W;% time-bandwidth product
k = floor(2*(T*W)-1);% number of tapers

params2              = [];
params2.Fs           = Fs;
params2.pad          = 4;
params2.fpass        = [20 220];
params2.tapers       = [TW k];
params2.trialave     = 1;

%% pre-calculate filter coefficients using finite impulse response filter (linear phase)
bpf = [2.5 12];
[IEDf1] = fir1(3*floor(Fs/bpf(1)),[ bpf ]./(Fs/2),'bandpass');

bpf = [60 120];
[IEDf2] = fir1(3*floor(Fs/bpf(1)),[ bpf ]./(Fs/2),'bandpass');

bpf = [2.5];
[HPf1] = fir1(3*floor(Fs/bpf(1)),[ bpf ]./(Fs/2),'high');

%%
[rpath] = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',pID,'/',expMode,'/'];

[savePath] = [rpath(1:regexp(rpath,expMode)-1),'res/'];
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

%% load spk and LFP data
for seshIt = 1:length( sesh)
    
    [p2d]       = [rpath,sesh{seshIt},'/'];
    [spkDat]    = load([p2d,pID,'_',expMode,'_',sesh{seshIt},'_spkDataStimLockedSegmented.mat']);
    [lfpDat]    = load([p2d,pID,'_',expMode,'_',sesh{seshIt},'_lfpDataStimLockedSegmenteddownsampled.mat']);
    
    %% time range of interest and fS
    if isfield(lfpDat,'trlSel')
         [tIx] = find(lfpDat.dsTrlTime >=-1 & lfpDat.dsTrlTime <=4);
         [Fs] = lfpDat.dsFs;
    else% backward compatibility
        [tIx] = find(lfpDat.trlTime >=-1 & lfpDat.trlTime <=4);
        [Fs] = lfpDat.Fs;
    end;
 
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
    params4.fpass        = [20 220];
    params4.tapers       = [TW k];
    params4.trialave     = 1;
    params4.err          = [2 0.001];%jacknife error
    
    paramsTrl            = params3;
    paramsTrl.trialave   = 0;
    
    %% extract BF labels and ids
    ix = regexp(lfpDat.chanLab,'\d{1}');
    ix = [ix{:}]-1;
    [BFlab] = cell(1,length(ix));
    for it = 1:length(ix)
        BFlab(it) = { lfpDat.chanLab{it}(1:ix(it)) };
    end;
    [BFid] = unique(BFlab)';
    
    %% estimate median LFP & power spectrum
%     [AVG] = cell(1,length(BFid));
%     [ICs1] = cell(1,length(BFid));
%     [ICs2] = cell(1,length(BFid));
    [S1] = cell(1,length(BFid));
    [S2] = cell(1,length(BFid));
    [Serr1] = cell(1,length(BFid));
    [Serr2] = cell(1,length(BFid));
%     [weights1] = cell(1,length(BFid));
%     [icS1] = cell(length(BFid),length(BFid));
%     [icSerr1] = cell(length(BFid),length(BFid));
%     [weights2] = cell(1,length(BFid));
%     [icS2] = cell(length(BFid),length(BFid));
%     [icSerr2] = cell(length(BFid),length(BFid));
    
    %% recruit workers for parallel computing
    if isempty(gcp('nocreate'))
        parpool(36,'SpmdEnabled', false)
    end;
    
    %% trial indexes: hits vs misses
    [ntrl] = size(lfpDat.LFPseg{1},2);
    
    [missIdx2] = find(ismember(sort([lfpDat.missIdx;lfpDat.hitIdx]),lfpDat.missIdx));
    [hitIdx2] = find(ismember(sort([lfpDat.missIdx;lfpDat.hitIdx]),lfpDat.hitIdx));
    
    xx = zeros(ntrl,1);
    xx(hitIdx2) = 1;% keep only hits
    
    %%
    [ selTrl ]   = cell(length(BFlab),1);
    [ seTrlAVG ] = cell(length(BFid),1);
    for it = 1:length(BFid)
        
        fprintf([num2str(it),'/',num2str(length(BFid)),'\n']);
        
        %%
        ix = find( strcmp(BFlab,BFid{it}) );
        [nsmp,ntrl] = size( lfpDat.LFPseg{ix(1)}(:,xx==1));
        
        %% reject trials with dubious median LFP activity (>4SDs)
        for jt = 1:length(ix)
            
            lfpDat.LFPseg{ix(jt)} = lfpDat.LFPseg{ix(jt)}(:,xx==1);
            
            x = lfpDat.LFPseg{ix(jt)};
            
            filt1 = zeros(size(x));
            filt2 = zeros(size(x));
            parfor kt = 1:size(x,2)
                d = x(:,kt)';
                d1 = abs(hilbert(filtfilt(IEDf1,1,[fliplr(d) d fliplr(d)]))).^2;
                d2 = abs(hilbert(filtfilt(IEDf2,1,[fliplr(d) d fliplr(d)]))).^2;
                
                filt1(:,kt) = d1(length(d)+1:length(d)*2);
                filt2(:,kt) = d2(length(d)+1:length(d)*2);
            end;
            
            filt1 =  filt1(tIx,:);
            M = ones(size(filt1,1),1)*mean(filt1,1);
            SD = ones(size(filt1,1),1)*std(filt1,0,1);
            z1 = (filt1-M)./SD;
            z1 = max(abs(z1),[],1);
            
            filt2 = filt2(tIx,:);
            M = ones(size(filt2,1),1)*mean(filt2,1);
            SD = ones(size(filt2,1),1)*std(filt2,0,1);
            z2 = (filt2-M)./SD;
            z2 = max(abs(z2),[],1);
            
            x =  x(tIx,:);
            M = ones(size(x,1),1)*mean(x,1);
            SD = ones(size(x,1),1)*std(x,0,1);
            z3 = (x-M)./SD;
            z3 = max(abs(z3),[],1);
            
            lfpDat.LFPseg{ix(jt)} = lfpDat.LFPseg{ix(jt)}(:,z1<8 & z2<8 & z3<4);% reject ARTEFACTS
            %lfpDat.LFPseg{ix(jt)} = lfpDat.LFPseg{ix(jt)}(:,z3<4);% reject ARTEFACTS
            
            selTrl{ix(jt)} = find(z1<8 & z2<8 & z3<4);
            %selTrl{ix(jt)} = find(z3<4);
        end;
        clear filt* M SD z*;
        
        %%
        dx = min([selTrl{ix}]):max([selTrl{ix}]);
        n= zeros(length(ix),length(dx));
        parfor jt = 1:length(ix)
            [n(jt,:)] = hist(selTrl{ix(jt)},dx);
        end;
        [trlSel] = dx(find(sum(n,1)==length(ix)));
        
        [ seTrlAVG{it} ] = trlSel;
        
        %% demean, detrend & high-pass filter individual trials of LFP data
        for jt = 1:length(ix)
            x = lfpDat.LFPseg{ix(jt)};
            parfor kt = 1:size(x,2)
                d =locdetrend( x(:,kt),Fs,[1 .25]);
                d = d -mean(d);% mean over time
                %f = filtfilt(HPf1,1,[fliplr(d') d' fliplr(d')]);
                %f = f(length(x(:,kt))+1:length(x(:,kt))*2);
                x(:,kt) = d;
            end;
            lfpDat.LFPseg{ix(jt)} = x;
        end;
        clear x;
        
%         %% bandpass filter and concatenate
%         [concat1] = zeros(length(ix),length(tIx)*length(trlSel));
%         [concat2] = zeros(length(ix),length(tIx)*length(trlSel));
%         ixi = 1:length(tIx);
%         for kt = 1:length( trlSel )
%             for jt = 1:length(ix)
%                 x = lfpDat.LFPseg{ix(jt)}(:,ismember( selTrl{ix(jt)}, trlSel(kt)));
%                 ft = filtfilt(b5,1,[fliplr(x) x fliplr(x)]);% filter <30Hz
%                 ft = ft(nsmp+1:2*nsmp);
%                 concat1(jt,ixi) = ft(tIx)-mean(ft(tIx));% demean
%                 ft = filtfilt(b4,1,[fliplr(x) x fliplr(x)]);% filter >30Hz
%                 ft = ft(nsmp+1:2*nsmp);
%                 concat2(jt,ixi) = ft(tIx)-mean(ft(tIx));% demean
%             end;
%             ixi = ixi+length(tIx);
%         end;
%         
%         %     %% PCA
%         %     [u,s,v] = svd(concat1,'econ');
%         %     sel = cumsum(diag(s)./sum(diag(s)))<=0.98;
%         %     [rec1] = u(:,sel)*s(sel,sel)*v(:,sel)';
%         %
%         %     [u,s,v] = svd(concat2,'econ');
%         %     sel = cumsum(diag(s)./sum(diag(s)))<=0.98;
%         %     [rec2] = u(:,sel)*s(sel,sel)*v(:,sel)';
%         
%         %% ica decomposition
%         %nComp = rank(concat1*concat1');
%         %[weights1{it},sphere1] = runica(concat1,'pca',nComp);
%         [weights1{it},sphere1] = runica(concat1);
%         
%         %nComp = rank(concat2*concat2');
%         %[weights2{it},sphere2] = runica(concat2,'pca',nComp);
%         [weights2{it},sphere2] = runica(concat2);
%         
%         datamean1 = mean(concat1')';
%         datamean2 = mean(concat2')';
%         [ICs1{it}] = icaact(concat1,weights1{it}*sphere1,datamean1);%compute timeseries for ICs
%         [ICs2{it}] = icaact(concat2,weights2{it}*sphere2,datamean2);%compute timeseries for ICs
%         clear datamean* sphere* recon* concat*;
%         
%         %% compute correlation across channels
%         [x] = zeros(length(ix),nsmp,length( trlSel ));
%         for jt = 1:length( ix )
%             y = zeros(length(trlSel),nsmp);
%             parfor kt = 1:length( trlSel )
%                 d = lfpDat.LFPseg{ix(jt)}(:,ismember( selTrl{ix(jt)}, trlSel(kt)));%select HITS
%                 d = filtfilt(b5,1,[fliplr(d) d fliplr(d)]);
%                 d = d(nsmp+1:2*nsmp);
%                 y(kt,:)= d;
%             end;
%             x(jt,:,:) = y';
%         end;
%         clear d y;
%         
%         rxy = zeros(size(x,1),size(x,1));
%         for jt = 1:size(x,1)
%             for kt = 1:size(x,1)
%                 rxy(jt,kt) = mean(diag(corr(squeeze(x(jt,:,:)),squeeze(x(kt,:,:)))));
%                 if jt==kt
%                     rxy(jt,kt) = 0;
%                 end;
%             end;
%         end;clear x;
%         rxy =(rxy>.7);
%         ix = ix(find(sum(rxy,2) >= median(sum(rxy,2))));
%         clear rxy x;
%         
%         %% compute median LFP across channels that are sufficiently highly correlated (rxy>.7)
%         x = zeros(length(ix),nsmp,length(trlSel));
%         for kt = 1:length( trlSel )
%             for jt = 1:length( ix )
%                 x(jt,:,kt) = lfpDat.LFPseg{ix(jt)}(:,ismember( selTrl{ix(jt)}, trlSel(kt)));%select HITS
%             end;
%         end;
%         x = squeeze(median(x,1));
%         AVG{it} = x;
%         clear x;
%         
%         %% demean median LFP
%         AVG{it} = AVG{it}-ones(size(AVG{it},1),1)*mean(AVG{it},1);
%         
%         %% compute power spectrum
%         [S1{it},f1,Serr1{it}] = mtspectrumc( gradient(AVG{it}(tIx,:)')', params3 );
%         [S2{it},f2,Serr2{it}] = mtspectrumc( gradient(AVG{it}(tIx,:)')', params4 );
%         
%         %%
%         [dum1] = reshape(ICs1{it},[size(ICs1{it},1) length(tIx) length(trlSel)]);
%         
%         %res = 1/((2^(nextpow2(length(tIx))+params3.pad))/Fs);
%         %[nfft] = length(params3.fpass(1)+res/1.55:res:params3.fpass(2));
%         sx = zeros(length(f1),size( ICs1{it},1));
%         sxe = zeros(2,length(f1),size( ICs1{it},1));
%         nICs = size( ICs1{it},1);
%         parfor jt = 1:size( ICs1{it},1)
%             fprintf([num2str(jt),'/',num2str(nICs)]);
%             X = squeeze(dum1(jt,:,:));
%             X = X-ones(size(X,1),1)*mean(X,1);
%             [sx(:,jt),~,sxe(:,:,jt)] = mtspectrumc( gradient( X' )' , params3 );
%             fprintf('\n');
%         end;
%         icS1{it}    = sx;
%         icSerr1{it} = sxe;
%         clear dum* sx*;
%         
%         [dum2] = reshape(ICs2{it},[size(ICs2{it},1) length(tIx) length(trlSel)]);
%         
%         %res = 1/((2^(nextpow2(length(tIx))+params4.pad))/Fs);
%         %[nfft] = length(params4.fpass(1)+res/1.55:res:params4.fpass(2));
%         sx = zeros(length(f2),size( ICs1{it},1));
%         sxe = zeros(2,length(f2),size( ICs1{it},1));
%         parfor jt = 1:size( ICs2{it},1)
%             fprintf([num2str(jt),'/',num2str( nICs )]);
%             X = squeeze(dum2(jt,:,:));
%             X = X-ones(size(X,1),1)*mean(X,1);
%             [sx(:,jt),~,sxe(:,:,jt)] = mtspectrumc( gradient( X' )' , params4 );
%             fprintf('\n');
%         end;
%         icS2{it} = sx;
%         icSerr2{it} = sxe;
%         clear dum* sx*;
        
    end;
    
    %% save data
    saveName = [pID,'_',expMode,'_',sesh{seshIt},'_preprocLFP.mat'];
    save([savePath,saveName],'lfpDat','selTrl');
%     
%     saveName = [pID,'_',expMode,'_',sesh{seshIt},'_medianANDicaLFP.mat'];
%     save([savePath,saveName],'AVG','ICs1','weights1','ICs2','weights2','seTrlAVG');
%     clear ICs;
%     
%     saveName = [pID,'_',expMode,'_',sesh{seshIt},'_medianLFP_powerSpectra.mat'];
%     save([savePath,saveName],'S1','S2','f1','f2','Serr1','Serr2');
%     
%     saveName = [pID,'_',expMode,'_',sesh{seshIt},'_LFPICs_powerSpectra.mat'];
%     save([savePath,saveName],'icS1','icS2','f1','f2','icSerr1','icSerr2');
%     
%     %% visualize spectra
%     h = figure;
%     a  = zeros(length(AVG),1);
%     lm = zeros(length(AVG),4);
%     for it = 1:length( AVG )
%         subplot(length(AVG)/2,2,it);
%         a(it) = gca;
%         jbfill(f1,Serr1{it}(1,:),Serr1{it}(2,:),[0 0 .9],[0 0 .9],1,.5);
%         jbfill(f2,Serr2{it}(1,:),Serr2{it}(2,:),[.9 0 0],[.9 0 0],1,.5);
%         hold on;
%         plot(f1,S1{it},'b');
%         plot(f2,S2{it},'r');
%         lm(it,:) = [min(S1{it}) max(S1{it}) min(S2{it}) max(S2{it})];
%         title(BFid{it});
%         xlabel('Frequency (Hz)');
%         ylabel('Power (a.u.)');
%     end;
%     set(a,'YLim',[min(min(lm)) max(max(lm))]);
%     
%     saveName = [pID,'_',expMode,'_',sesh{seshIt},'_medianLFP_powerSpectra.fig'];
%     saveas(h,[savePath,saveName]);
%     clear S1 S2 Serr1 Serr2;
%     
%     %% visualize IC-spectra
%     C = {'r','g','b','c','m','y','k',[.5 .5 .5]};
%     
%     for it = 1:length(weights1)
%         
%         figure;
%         subplot(221);
%         hold on;
%         for kt = 1:size(weights1{it},2)
%             plot(weights1{it}(:,kt),1:size(weights1{it},1),'Color',C{kt})
%         end;
%         plot([0 0],[0 size(weights1{it},2)+1],'k--')
%         axis tight;xlim([-1 1]);
%         title('Spatial Weights');
%         ylabel('mw#');
%         xlabel('(a.u.)');
%         
%         subplot(222);
%         for kt = 1:length(icS1)
%             %jbfill(f1,icSerr1{it}(1,:),icSerr1{it}(2,:),C{kt},C{kt},1,.5);
%             hold on;
%             plot(f1,icS1{it}(:,kt),'Color',C{kt});
%         end;
%         xlabel('Frequency (Hz)');
%         ylabel('Power (a.u.)');
%         xlim([0 20]);
%         
%         subplot(223);
%         hold on;
%         for kt = 1:size(weights2{it},2)
%             plot(weights2{it}(:,kt),1:size(weights2{it},1),'Color',C{kt})
%         end;
%         plot([0 0],[0 size(weights2{it},2)+1],'k--')
%         axis tight;xlim([-1 1]);
%         title('Spatial Weights');
%         ylabel('mw#');
%         xlabel('(a.u.)');
%         
%         subplot(224);
%         for kt = 1:length(icS2)
%             %jbfill(f2,icSerr2{it}(1,:),icSerr2{it}(2,:),C{kt},C{kt},1,.5);
%             hold on;
%             plot(f2,icS2{it}(:,kt),'Color',C{kt});
%         end;
%         xlabel('Frequency (Hz)');
%         ylabel('Power (a.u.)');
%         xlim([40 220]);
%         
%         saveName = [pID,'_',expMode,'_',sesh{seshIt},'_',BFid{it},'_LFPICs_powerSpectra.fig'];
%         saveas(h,[savePath,saveName]);
%     end;
%     clear weights* icS* icSerr*;
%     
%     %% compute coherence spectra lower freqz
%     T = floor(length(tIx)/Fs);
%     W = 3;
%     TW = T*W;
%     k = 2*TW-1;
%     
%     paramsC                     = [];
%     paramsC.Fs                  = Fs;
%     paramsC.tapers              = [TW k];
%     paramsC.pad                 = 2;
%     paramsC.fpass               = params3.fpass;
%     paramsC.err                 = [2 0.001];
%     paramsC.trialave            = 1;
%     
%     npairs = length(BFid)*(length(BFid)-1)/2;
%     p = zeros(npairs,2);
%     cnt = 0;
%     for it = 1:length(BFid)
%         for jt = it+1:length(BFid)
%             cnt = cnt+1;
%             p(cnt,:) = [it jt];
%         end;
%     end;
%     
%     [C1] = cell( 1, npairs );
%     [Cerr1] = cell( 1, npairs );
%     np = size(p,1);
%     parfor it = 1:size(p,1)
%         fprintf([num2str(it),'/',num2str(np)]);
%         seTrlAVG2 = seTrlAVG;
%         X = AVG{p(it,1)}(tIx, ismember(seTrlAVG2{p(it,1)},intersect(seTrlAVG2{p(it,1)},seTrlAVG2{p(it,2)})) );
%         Y = AVG{p(it,2)}(tIx, ismember(seTrlAVG2{p(it,2)},intersect(seTrlAVG2{p(it,1)},seTrlAVG2{p(it,2)})) );
%         [C1{it},~,~,~,~,~,~,~,Cerr1{it}] = coherencyc( X, Y, paramsC );
%         fprintf('\n');
%     end;
%     
%     saveName = [pID,'_',expMode,'_',sesh{seshIt},'_medianLFP_coherenceSpectra.mat'];
%     save([savePath,saveName],'C1','f1','Cerr1');
%     
%     %% compute coherence spectra higher freqz
%     T = floor(length(tIx)/Fs);
%     W = 8;
%     TW = T*W;
%     k = 2*TW-1;
%     
%     paramsC                     = [];
%     paramsC.Fs                  = Fs;
%     paramsC.tapers              = [TW k];
%     paramsC.pad                 = 2;
%     paramsC.fpass               = params4.fpass;
%     paramsC.err                 = [2 0.001];
%     paramsC.trialave            = 1;
%     
%     npairs = length(BFid)*(length(BFid)-1)/2;
%     p = zeros(npairs,2);
%     cnt = 0;
%     for it = 1:length(BFid)
%         for jt = it+1:length(BFid)
%             cnt = cnt+1;
%             p(cnt,:) = [it jt];
%         end;
%     end;
%     
%     [C2] = cell( 1, npairs );
%     [Cerr2] = cell( 1, npairs );
%     np = size(p,1);
%     parfor it = 1:size(p,1)
%         fprintf([num2str(it),'/',num2str(np)]);
%         seTrlAVG2 = seTrlAVG;
%         X = AVG{p(it,1)}(tIx, ismember(seTrlAVG2{p(it,1)},intersect(seTrlAVG2{p(it,1)},seTrlAVG2{p(it,2)})) );
%         Y = AVG{p(it,2)}(tIx, ismember(seTrlAVG2{p(it,2)},intersect(seTrlAVG2{p(it,1)},seTrlAVG2{p(it,2)})) );
%         [C2{it},~,~,~,~,~,~,~,Cerr2{it}] = coherencyc( X, Y, paramsC );
%         fprintf('\n');
%     end;
%     clear AVG;
%     
%     %% save coherence data
%     saveName = [pID,'_',expMode,'_',sesh{seshIt},'_medianLFP_coherenceSpectra.mat'];
%     save([savePath,saveName],'C1','f1','Cerr1','C2','f2','Cerr2');
%     
%     %%
%     pid = unique(p(:,1));
%     cnt2 = 0;
%     h = figure;
%     for it = 1:length(pid)
%         ix = find(p(:,1)==pid(it));
%         cnt = ((it-1)*length(pid))+(length(pid)-length(ix));
%         for jt = 1:length(ix)
%             cnt = cnt+1;
%             cnt2 = cnt2+1;
%             subplot(length(pid),length(pid),cnt)
%             jbfill(f1,Cerr1{cnt2}(1,:),Cerr1{cnt2}(2,:),[0 0 .9],[0 0 .9],1,.5);
%             jbfill(f2,Cerr2{cnt2}(1,:),Cerr2{cnt2}(2,:),[.9 0 0],[.9 0 0],1,.5);
%             hold on;
%             plot(f1,C1{cnt2},'y');
%             plot(f2,C2{cnt2},'y');
%             axis tight;
%             title(BFid{p(ix(jt),2)});
%             ylabel(BFid{p(ix(jt),1)});
%         end;
%     end;
%     
%     saveName = [pID,'_',expMode,'_',sesh{seshIt},'_medianLFP_coherenceSpectra.fig'];
%     saveas(h,[savePath,saveName]);
%     clear f1 f2 C* Cerr*;
    
    %%            
    if ~isempty(gcp('nocreate'))
        delete(gcp);
    end;       
    
    %%
    tic;
    for it = 1:length(BFlab) % loop over the BFs
        
        fprintf([num2str(it),'/',num2str(length(BFlab))]);
        
        %%
        if length( selTrl{it} ) >= 20
            
            %%
            [nsmp,ntrl] = size(lfpDat.LFPseg{it});% get dims of input
            
            %% compute power spectra
            [S1,f1] = mtspectrumc( gradient(lfpDat.LFPseg{it}(tIx,:)')', paramsTrl);
            [~,mIx] = max(mean(S1,2));
            pfL = f1(mIx);
            powPFL = mean(S1(f1 >= f1(mIx)-2 & f1 <=f1(mIx)+1,:),2);
            [~,pvalPowL] = ttest(powPFL);
            
            [S1,f1,Serr1] = mtspectrumc( gradient(lfpDat.LFPseg{it}(tIx,:)')', params3);
            
            [S2,f2,Serr2] = mtspectrumc( gradient(lfpDat.LFPseg{it}(tIx,:)')', params4);
            [~,mIx] = max(S1);
            pfH = f2(mIx);
            
            [S3,t3,f3] = mtspecgramc( gradient(lfpDat.LFPseg{it}')', movingwin1, params1);
            t3 = t3-5;
            
	    if isfield(lfpDat,'dsTrlTime')
		S3 = S3(t3>=min(lfpDat.dsTrlTime(tIx)) & t3<=max(lfpDat.dsTrlTime(tIx)),:,:);
            	t3 = t3(t3>=min(lfpDat.dsTrlTime(tIx)) & t3<=max(lfpDat.dsTrlTime(tIx)));
	    else
            	S3 = S3(t3>=min(lfpDat.trlTime(tIx)) & t3<=max(lfpDat.trlTime(tIx)),:,:);
            	t3 = t3(t3>=min(lfpDat.trlTime(tIx)) & t3<=max(lfpDat.trlTime(tIx)));
            end;

            [S4,t4,f4] = mtspecgramc( gradient(lfpDat.LFPseg{it}')', movingwin2, params2);
            t4 = t4-5;

	    if isfield(lfpDat,'dsTrlTime')
		S4 = S4(t4>=min(lfpDat.dsTrlTime(tIx)) & t4<=max(lfpDat.dsTrlTime(tIx)),:);
            	t4 = t4(t4>=min(lfpDat.dsTrlTime(tIx)) & t4<=max(lfpDat.dsTrlTime(tIx)));
	    else
            	S4 = S4(t4>=min(lfpDat.trlTime(tIx)) & t4<=max(lfpDat.trlTime(tIx)),:);
            	t4 = t4(t4>=min(lfpDat.trlTime(tIx)) & t4<=max(lfpDat.trlTime(tIx)));
            end;

            saveName = [pID,'_',expMode,'_',sesh{seshIt},'_',lfpDat.chanLab{it},'_spectralData.mat'];
            save([savePath,saveName],'S1','f1','Serr1','S2','f2','Serr2','S3','t3','f3','S4','t4','f4','pfL','pfH','pvalPowL','powPFL');
            
            %%            
            h = [];
            h(1) = figure;
            hold on;
            jbfill(f1,(Serr1(1,:)),(Serr1(2,:)),[0 0 .9],[0 0 .9],1,.5);
            jbfill(f2,(Serr2(1,:)),(Serr2(2,:)),[.9 0 0],[.9 0 0],1,.5);
            hold on;
            plot(f1,(S1),'y');
            plot(f2,(S2),'y');
            title(lfpDat.chanLab{it});
            xlabel('Frequency (Hz)');
            ylabel('Power (a.u.)');
            
            h(2) = figure;
            subplot(6,1,1:4);
            hold on;
            imagesc(t4.*1e3,f4,20*log10(S4)');
            plot([0 0],[min(f4) max(f4)],'w-');
            plot([2e3 2e3],[min(f4) max(f4)],'w-');
            axis xy;axis tight;
            ylabel('Frequency (Hz)');
            xlabel('Time (ms)');
            subplot(6,1,5:6);
            hold on;
            imagesc(t3.*1e3,f3,20*log10(S3)');
            plot([0 0],[min(f3) max(f3)],'w-');
            plot([2e3 2e3],[min(f3) max(f3)],'w-');
            axis xy;axis tight;
            ylabel('Frequency (Hz)');
            xlabel('Time (ms)');
            
            for kt = 1:length(h)
                saveName = [pID,'_',expMode,'_',sesh{seshIt},'_',lfpDat.chanLab{it},'_spectralData',num2str(kt),'.fig'];
                saveas(h(kt),[savePath,saveName],'fig');
            end;
            clear S1 S2 f1 f2 Serr1 Serr2 S3 t3 f3 S4 t4 f4;
            
        end;
        
        %%
        close all;
        fprintf('\n');
        
    end;
    toc;
    
end;
