%%
addpath(genpath('/home/rouxf/tbx/chronux_2_11/'));
addpath(genpath('/home/rouxf/prj/Bham/code/mcode/utils/'));

%%
T = 1;
W = 3;
TW = T*W;
k = 2*TW-1;

params1              = [];
params1.Fs           = 1e3;
params1.pad          = 0;
params1.fpass        = [1 30];
params1.tapers       = [TW k];
params1.trialave     = 1;
params1.err          = [1 0.001];

T = 1;
W = 5;
TW = T*W;
k = 2*TW-1;

params2              = [];
params2.Fs           = 1e3;
params2.pad          = 0;
params2.fpass        = [20 170];
params2.tapers       = [TW k];
params2.trialave     = 1;
params2.err          = [1 0.001];

%%
[pID]     = 'P07';
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
SEG = cell(1,n);
for it = 1:length(sesh)
    
    [path2d] = [rpath,pID,'/',expMode,'/',sesh{it},'/'];
    [fn] = dir([path2d,pID,'_',expMode,'_',sesh{it},'_lfpDataStimLockedContinuousDownsampled.mat']);
    
    if ~isempty( fn )
        load([path2d,fn.name]);
    end;
    
    %%
    lfpAVG = cell(1,length(chanLab));
    ix = 1:8;
    for jt = 1:length( chanLab )
        
        x = zeros(size(LFPsig{1}));
        for kt = 1:length( ix )
            x = x + LFPsig{ix(kt)};
        end;
        x = x./kt;
        
        [lfpAVG{jt}] = x;
        clear x;
        ix = ix+8;
    end;
    
    %%
    n = Fs;
    m = length(lfpAVG{1});
    p = floor(m/n);
    
    seg = cell(1,length( chanLab ));
    for kt = 1:length( chanLab )
        fprintf([num2str(kt),'/',num2str(length(chanLab))]);
        ix = 1:Fs;
        for jt = 1:p
            seg{kt}(:,jt) = lfpAVG{kt}(ix);
            ix = ix+Fs;
        end;
        SEG{kt} = [SEG{kt}';seg{kt}']';
        clear seg;
        fprintf('\n');
    end;
    
end;
SEG = SEG(1:length( chanLab) );
clear lfpAVG;

%%
dy =  cell( 1, length( chanLab ) );
for jt = 1:length( chanLab )
    
    [n,m] = size(SEG{jt});
    x = SEG{jt};
    m2 = ones(size(x,1),1)*(max(x,[],1));
    m1 = ones(size(x,1),1)*(min(x,[],1));
    mx = ones(n,1)*mean(x,1);
    sx = ones(n,1)*std(x,0,1);
    z = (x-mx)./sx;
    s = max(abs(z),[],1);
    %x = (x-m1)./(m2-m1);
    SEG{jt} = z(:,s<3);% removes all trials with SD>4
    
    [dy{jt}] = gradient( SEG{jt}' )';  % dx = horizontal, dy = vertical
    %dy{jt} = SEG{jt};
end;
clear x m1 m2 mx sx z s SEG;
    
%%
S1 = zeros(length(dy),29);
S2 = zeros(length(dy),154);
SE1 = zeros(length(dy),2,29);
SE2 = zeros(length(dy),2,154);
for jt = 1:length( dy )
    fprintf([num2str(jt),'/',num2str( length( dy ) )]);
    [S1(jt,:),fx1, SE1(jt,:,:)] = mtspectrumc( dy{jt}, params1);    
    [S2(jt,:),fx2, SE2(jt,:,:)] = mtspectrumc( dy{jt}, params2);
    fprintf('\n');
end;
clear dy;
  
% %%
% c = zeros(size(S1));
% for it = 1:size(S1,1)
%     
%     X = log10(fx1)';
%     Y = log10(S1(it,:))';
%     
%     [b] = regress(Y,[ones(size(X)) X]);
%     
%     yp = [ones(size(X)) X]*b;
%     
%     d = Y-yp;
%     c(it,:) = 10.^d;
%     
% end;
% S1 = c;
% 
% c = zeros(size(S2));
% for it = 1:size(S2,1)
%     
%     X = log10(fx2)';
%     Y = log10(S2(it,:))';
%     
%     [b] = regress(Y,[ones(size(X)) X]);
%     
%     yp = [ones(size(X)) X]*b;
%     
%     d = Y-yp;
%     c(it,:) = 10.^d;
%     
% end;
% S2 = c;

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
subplot(221);
hold on;
for it = 1:length(ixL)
    plot(fx1,S1(sIx(ixL(it)),:),'r','LineWidth',3);
    jbfill(fx1,squeeze(SE1(sIx(ixL(it)),1,:))',squeeze(SE1(sIx(ixL(it)),2,:))','r','r');
    hold on;
end;
axis tight;
subplot(223);
hold on;
for it = 1:length( ixR )
    plot(fx1,squeeze(median(S1(sIx(ixR(it)),:),1)),'b','LineWidth',3);
    jbfill(fx1,squeeze(median(SE1(sIx(ixR(it)),1,:),1))',squeeze(median(SE1(sIx(ixR(it)),2,:),1))','b','b'); 
    hold on;
end;
axis tight;

subplot(222);
hold on;
for it = 1:length( ixL )
    plot(fx2,S2(sIx(ixL(it)),:),'r','LineWidth',3);
    jbfill(fx2,squeeze(SE2(sIx(ixL(it)),1,:))',squeeze(SE2(sIx(ixL(it)),2,:))','r','r');
    hold on;
end;
axis tight;
subplot(224);
hold on;
for it = 1:length( ixR )
    plot(fx2,S2(sIx(ixR(it)),:),'b','LineWidth',3);
    jbfill(fx2,squeeze(SE2(sIx(ixR(it)),1,:))',squeeze(SE2(sIx(ixR(it)),2,:))','b','b');
    hold on;
end;
axis tight;