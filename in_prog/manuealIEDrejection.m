%%
T = size(LFPavg{1},1)/Fs;
W = 1/T;
TW = round(T*W);
k = round(2*TW-1);
params                  = [];
params.Fs               = Fs;
params.pad              = 3;
params.fpass            = [40 60];
params.trialave         = 0;
params.tapers           = [TW k];

for jt = 1:length(LFPavg)
    fprintf([num2str(jt),'/',num2str(length(LFPavg))]);
    for it = 1:size(LFPavg{jt},2)
        
        x = LFPavg{jt}(:,it);
        [S,f] = mtspectrumc( x-mean(x), params );
        
        S = 20*log10(S);
        S = (S-mean(S))./std(S);
        [~,ix] = max(S);
        
        %[b,a] = butter(4,[f(ix)-.5 f(ix)+.5]./(Fs/2),'stop');% apply band-stop for LFP
        %[dum] = filtfilt(b,a,LFPavg{jt}(:,it)');
        
        [dum] = CleanLineNoise(LFPavg{jt}(:,it)','Fs',Fs,'noiseFreq',f(ix),'windowSize',T);
        
%         x1 = [];
%         x2 = [];
%         t = [1:length(x)]./Fs;
%         x1 = sin(2*pi*f(ix).*t)';
%         x2 = cos(2*pi*f(ix).*t)';
%         
%         X = [x1 x2 x1+x2 x1.*x2];
%         
%         Y =x;
%         Y = Y-mean(Y);
%         
%         b = regress(Y,[ones(size(X,1),1) X]);
%         
%         yp = [ones(size(X,1),1) X]*b;
        
        %LFPavg{jt}(:,it) = x-yp;
        LFPavg{jt}(:,it) = dum;
        
        %     x = dum;
        %     [S2,f] = mtspectrumc( (x-mean(x))', params );
        %     S2 = 20*log10(S2);
        %     S2 = (S2-mean(S2))./std(S2);
        %
        %     figure;
        %     subplot(121);
        %     plot(f,S);hold on;
        %     plot(f,S2,'r');
        %
        %     dum=rmlinesc(x,params,0.05,'n',f(ix));
        %     x = dum;
        %     [S2,f] = mtspectrumc( (x-mean(x))', params );
        %     S2 = 20*log10(S2);
        %     S2 = (S2-mean(S2))./std(S2);
        %
        %     subplot(122);
        %     plot(f,S);hold on;
        %     plot(f,S2,'r');
    end;
    fprintf(['\n']);
end;

%%
dum = [];
dum.label = chanLab';

ntrl = size(LFPavg{1},2);
for jt = 1:ntrl
    for it = 1:length( LFPavg )
        x = LFPavg{it}(:,jt);
        dum.trial{jt}(it,:) = x;
        dum.time{jt}        = trlTime;
        
    end;
end;

cfg                     = [];
cfg.demean              = 'yes';
cfg.detrend             = 'yes';
cfg.lpfilter            = 'yes';
cfg.lpfreq              = 100;

[dum] = ft_preprocessing( cfg, dum );

% for it = 1:ntrl
%         fprintf([num2str(it),'/',num2str(ntrl)]);
%         X = dum.trial{it};
%         M  = mean(X,2);
%         SD = std(X,0,2);
%         M=M*ones(1,size(X,2));
%         SD=SD*ones(1,size(X,2));
%         
%         dum.trial{it} = (X-M)./SD;
%         fprintf('\n');
% end;

cfg                     = [];
cfg.latency             = [-.5 4];

[dum] = ft_selectdata( cfg , dum );

cfg                     = [];
cfg.viewmode            = 'vertical';
%cfg.channel             = chanLab(5);

[cfg] = ft_databrowser( cfg, dum );

%%
nsamp = length(dum.trial{1});
ix = 1:nsamp;
trl = [];
for it = 1:length(dum.trial)
    trl(it,:) = ix+nsamp*(it-1);
end;

trl2 = cfg.artfctdef.visual.artifact;
sel = [];
sel2 = [];
for jt = 1:size(trl2,1)
    x = trl2(jt,1):trl2(jt,2);    
    for it = 1:size(trl,1)        
        if any(ismember(trl(it,:),x))
            sel = [sel it];
            sel2 = [sel2 jt];
        end;
    end;
end;

[delIx] = unique(sel);

%%
badChans = {'midHippR','postHippL'};
chanIx = [];
for it = 1:length(badChans)
    chanIx = [chanIx find(strcmp(chanLab,badChans{it}))];
end;

%%
for it = 1:length(LFPavg)   
    LFPavg{it}(:,delIx) = [];    
end;

LFPavg(chanIx) = [];

%%
movingwin1 = [1 0.001];
T = movingwin1(1);
W = 4;
TW = T*W;

params1                  =[];
params1.pad              = 0;
params1.Fs               = 1e3;
params1.fpass            = [0 30];
params1.tapers           = [TW 2*TW-1];
params1.trialave         = 1;

movingwin2 = [0.25 0.001];
T = movingwin2(1);
W = 10;
TW = T*W;

params2                  =[];
params2.pad              = 0;
params2.Fs               = 1e3;
params2.fpass            = [30 150];
params2.tapers           = [TW 2*TW-1];
params2.trialave         = 1;

Sl = cell(1,length(LFPavg));
Sh = cell(1,length(LFPavg));
for it = 1:length(LFPavg)
    fprintf([num2str(it),'/',num2str(length(LFPavg))]);
    
    dum = LFPavg{it};
    %M = repmat(mean(dum,1),[size(dum,1) 1]);
    %dum = dum-M;
    
    [Sl{it},tl,fl] = mtspecgramc( gradient(dum')', movingwin1, params1 );
    [Sh{it},th,fh] = mtspecgramc( gradient(dum')', movingwin2, params2 );
    
    fprintf('\n');
end;

%%
c = 0;
ix = [];
for it = 1:length( LFPavg )
    for jt = it+1:length( LFPavg )
      c = c+1;
      ix(c,:) = [it jt];
    end,
end;


T = length(LFPavg{1})/Fs;
W = 2;
TW = round(T*W);
k = 2*TW-1;

params                  = [];
params.Fs               = Fs;
params.pad              = 2;
params.fpass            = [0 30];
params.tapers           = [TW k]; 
params.trialave         = 1;

C = cell(size(ix,1),1);
for it = 1:size(ix,1)
    [C{it},~,~,~,~,fx] = coherencyc( LFPavg{ix(it,1)} , LFPavg{ix(it,2)}, params );
end;

T = length(LFPavg{1})/Fs;
W = 5;
TW = round(T*W);
k = 2*TW-1;

params2                  =[];
params2.pad              = 0;
params2.Fs               = 1e3;
params2.fpass            = [30 150];
params2.tapers           = [TW 2*TW-1];
params2.trialave         = 1;

C2 = cell(size(ix,1),1);
for it = 1:size(ix,1)
    [C2{it},~,~,~,~,fx2] = coherencyc( LFPavg{ix(it,1)} , LFPavg{ix(it,2)}, params2 );
end;

%%
ix = 1:6;
a = {};
ca1 = [];
ca2 = [];
for it = 1:length( LFPavg )
    figure;
    subplot(7,1,1:3);
    a{it} = gca;
    imagesc(th-5,fh,20*log10(Sh{it})');axis xy;
    %x = Sh{it}-ones(size(Sh{it},1),1)*mean(Sh{it}(th-5<0,:),1);
    %imagesc(th-5,fh,x');axis xy;
    xlim([-0.4 4]);
    ca1(it,:) = caxis;
    title(chanLab(ix(it)));
    subplot(7,1,5:6);
    a{it} = [a{it} gca];
    imagesc(tl-5,fl,20*log10(Sl{it})');axis xy;
    %x = Sl{it}-ones(size(Sl{it},1),1)*mean(Sl{it}(tl-5<50,:),1);
    %imagesc(tl-5,fl,x');axis xy;
    xlim([-0.4 4]);
    ca2(it,:) = caxis;
    subplot(7,1,7);
    hold on;
    x =LFPavg{it};
    plot(trlTime,mean(x,2),'b');
    axis tight;
    xlim([-0.4 4]);
    for jt = 1:length(a{it})
        hold(a{it}(jt),'on');
        plot(a{it}(jt),[0 0],[min(get(a{it}(jt),'YLim')) max(get(a{it}(jt),'YLim'))],'w');
        plot(a{it}(jt),[2 2],[min(get(a{it}(jt),'YLim')) max(get(a{it}(jt),'YLim'))],'w');
    end;
end;
   
% mL1 = [min(min(ca1)) max(max(ca1))];
% mL2 = [min(min(ca2)) max(max(ca2))];
% for it = 1:length(a)
%     caxis(a{it}(1),mL1);
%     caxis(a{it}(2),mL2);
% end;

%%
figure;
c = 0;k=0;
for it = 1:length(C)
    c=c+1;
    subplot(length(LFPavg)-1,length(LFPavg)-1,c);
    plot(fx,C{it});      
    axis tight;
    %ylim([0 1]);
    title([chanLab{ix(it,1)},'<->',chanLab{ix(it,2)}]);
    if ix(it,2)>=length(LFPavg)
        k = k+1;
        c = c+k;
    end;
end;

%%
figure;
c = 0;k=0;
for it = 1:length(C2)
    c=c+1;
    subplot(length(LFPavg)-1,length(LFPavg)-1,c);
    plot(fx2,C2{it});      
    axis tight;
    %ylim([0 1]);
    title([chanLab{ix(it,1)},'<->',chanLab{ix(it,2)}]);
    if ix(it,2)>=length(LFPavg)
        k = k+1;
        c = c+k;
    end;
end;