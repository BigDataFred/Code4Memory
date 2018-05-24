movingwin1 = [1 0.25];

T = movingwin1(1);
W = 2;
TW = T*W;
k = 2*TW-1;

params1              = [];
params1.Fs           = 1e3;
params1.pad          = 0;
params1.fpass        = [0 30];
params1.tapers       = [TW k];
params1.trialave     = 1;

movingwin2 = [0.25 0.067];

T = movingwin2(1);
W = 8;
TW = T*W;
k = 2*TW-1;

params2              = [];
params2.Fs           = 1e3;
params2.pad          = 0;
params2.fpass        = [30 170];
params2.tapers       = [TW k];
params2.trialave     = 1;

%%
[pID]     = 'P08';
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

%%
n = 8;
S1 = zeros(n,156,16);
S2 = zeros(n,158,46);
C1 = zeros(n*(n-1)/2,156,16);
C2 = zeros(n*(n-1)/2,158,46);
for it = 1:2%length(sesh)        
    
    [path2d] = [rpath,pID,'/',expMode,'/',sesh{it},'/'];
    [fn] = dir([path2d,pID,'_',expMode,'_',sesh{it},'_lfpDataStimLockedSegmentedAVGdownsampled.mat']);
    
    if ~isempty( fn )
        load([path2d,fn.name]);
    end;     
    
    %%
    p = [];
    c = 0;
    for jt = 1:length( chanLab )
        for kt = jt+1:length( chanLab )
            c = c+1;
            p(c,:) = [jt,kt];
        end;
    end;
    
    %%
    dy =  cell( 1, length( chanLab ) );
    for jt = 1:length( chanLab )
        
        [n,m] = size(LFPavg{jt});
        x = LFPavg{jt};
        mx = ones(n,1)*mean(x,1);
        sx = ones(n,1)*std(x,0,1);
        z = (x-mx)./sx;
        z = max(abs(z),[],1);
        LFPavg{jt} = LFPavg{jt}(:,z<4);% removes all trials with SD>4
        
        [dy{jt}] = gradient( LFPavg{jt}' )';  % dx = horizontal, dy = vertical
    end;
    %clear LFPavg;
    
    %%
    for jt = 1:length( dy )
        fprintf([num2str(jt),'/',num2str( length( dy ) )]);
        [s1,tx1,fx1] = mtspecgramc( dy{jt}, movingwin1, params1);
        tx1 = tx1-5;
        s1 = s1(tx1>=-1,:);
        tx1 = tx1(tx1>=-1);
        %s1 = 20*log10(s1);
        
        [s2,tx2,fx2] = mtspecgramc( dy{jt}, movingwin2, params2);
        tx2 = tx2-5;
        s2 = s2(tx2>=-1,:);
        tx2 = tx2(tx2>=-1);
        %s2 = 20*log10(s2);
        
        S1(jt,:,:) = squeeze(S1(jt,:,:)) + s1;
        S2(jt,:,:) = squeeze(S2(jt,:,:)) + s2;
        fprintf('\n');
    end;
    %clear dy;
    
    %%                
%     for jt = 1:size(p,1)
%         fprintf([num2str(jt),'/',num2str(size(p,1))]);
%         [c1,~,~,~,~,tc1,fc1] = cohgramc( dy{p(jt,1)}, dy{p(jt,2)}, movingwin1, params1 );
%         [c2,~,~,~,~,tc2,fc2] = cohgramc( dy{p(jt,1)}, dy{p(jt,2)}, movingwin2, params2 );
%         tc1 = tc1-5;
%         tc2 = tc2-5;
%         c1 = c1(tc1>=-1,:);
%         c2 = c2(tc2>=-1,:);
%         tc1 = tc1(tc1>=-1);
%         tc2 = tc2(tc2>=-1);
%         
%         C1(jt,:,:) = squeeze(C1(jt,:,:)) + c1;
%         C2(jt,:,:) = squeeze(C2(jt,:,:)) + c2;
%         fprintf('\n');
%     end;   
    
end;
S1 = S1./length(sesh);
S2 = S2./length(sesh);
C1 = C1./length(sesh);
C2 = C2./length(sesh);

n = length( chanLab);
S1 = S1(1:n,:,:);
S2 = S2(1:n,:,:);
C1 = C1(1:n*(n-1)/2,:,:);
C2 = C2(1:n*(n-1)/2,:,:);

%%
[hemLab] = cell( 1 , length( chanLab ) );
for it = 1:length( chanLab )
    hemLab(it) = { chanLab{it}(end) };
end;
[hemLab,sIx] = sort(hemLab);
[ixL] = find(strcmp(hemLab,'L'));
[ixR] = find(strcmp(hemLab,'R'));

%%
figure;
subplot(121);
plot(fx1,squeeze(median(mean(S1(sIx(ixL),:,:),2),1)),'r','LineWidth',3);
axis tight;
subplot(122);
plot(fx1,squeeze(median(mean(S1(sIx(ixR),:,:),2),1)),'LineWidth',3);
axis tight;

figure;
subplot(121);
plot(fx2,squeeze(median(mean(S2(sIx(ixL),:,:),2),1)),'r','LineWidth',3);
axis tight;
subplot(122);
plot(fx2,squeeze(median(mean(S2(sIx(ixR),:,:),2),1)),'LineWidth',3);
axis tight;

%%
figure;
for it = 1:size(S1,1)
    subplot(size(S1,1),1,it);
    imagesc(tx1,fx1,squeeze(S1(sIx(it),:,:))');
    axis xy;
    title(chanLab(sIx(it)));
end;

figure;
for it = 1:size(S1,1)
    subplot(size(S2,1),1,it);
    imagesc(tx2,fx2,squeeze(S2(sIx(it),:,:))');
    axis xy;
    title(chanLab(sIx(it)));
end;


%%
pID = unique(p(:,1));
c = 0;
figure;
for it = 1:length( pID )
    c = pID(it);
    ix = find(p(:,1) == pID(it));
    for jt = 1:length(ix)
        subplot(7,7,(jt+(it-1)*7)+pID(it)-1);
        hold on;
        plot(fc1,C1(ix(jt),:));
        plot(fc2,C2(ix(jt),:),'r');
        axis tight;ylim([0 1]);
        if jt == 1
            ylabel(chanLab{it});
        end;
        title(chanLab{p(ix(jt),2)});
    end;
    
end;



%%
ix = 1:8;
c = 0;
n = length(LFPsig)/8;
lfpAVG = cell(1,n);
while max(ix)<=length(LFPsig)
    c = c+1;
    fprintf([num2str(c),'/',num2str( n )]);
    x = zeros(1,length(LFPsig{1}));
    for jt = 1:length(ix)
        
        x = x + LFPsig{ ix(jt) };
        
    end;
    x = x./jt;
    ix = ix+8;
    
    lfpAVG{c} = x;
    fprintf('\n');
end;

%%
movingwin = [1 0.1];

T = movingwin(1);
W = 3;
TW = T*W;
k = 2*TW-1;

params              = [];
params.Fs           = 1e3;
params.pad          = 0;
params.fpass        = [0 30];
params.tapers       = [TW k];

S = cell( 1, length( lfpAVG ) );
for it = 1:length( lfpAVG )
    fprintf([num2str(it),'/',num2str( length( lfpAVG ) )]);
    [dy] = gradient(lfpAVG{it});
    [S{it},tx,fx] = mtspecgramc( dy', movingwin, params);
    fprintf('\n');
end;

%%
figure;
n = length( lfpAVG );
k = 0;
ix2 = 1:8;
for it = 1:length( lfpAVG )
    k =k+1;
    subplot(n,8,ix2(1:5));
    imagesc(tx,fx,20*log10(S{it})');
    axis xy;
    title(chanLab(it));
    xlabel('Time [s]');
    ylabel('Frequency [Hz]');
    subplot(n,8,ix2(end-1:end));
    plot(fx,mean(20*log10(S{it}),1),'LineWidth',3);
    axis tight;
    ix2 = ix2+8;
end;
colormap jet;