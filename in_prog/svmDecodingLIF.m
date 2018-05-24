%%
p2d = '/media/rouxf/rds-share/iEEG_DATA/MICRO/P07/fVSpEM/2017-04-30_17-35-52/';

load P07_fVSpEM_2017-04-30_17-35-52_spkDataStimLockedSegmented.mat;

%%
BFlab = cell(1,length(chanLab));
for it = 1:length( chanLab)
    
    x = chanLab{it};
    x(end) = [];
    BFlab(it) ={x};
end;
BFid = unique(BFlab)';

cIx = {[255 0 0]/255,[51 255 51]/255,[0 0 204]/255,[255 128 0]/255,[255 255 0]/255,[0 255 255]/255};

%%
cnt = 0;
wvf = {};
isi = {};
dt = 0:251;
selInfo = [];
wC = {};
for it = 1:length(sortedSpikes)
    
    [spk] = sortedSpikes{it};
    
    cID = unique(spk.assignedClusterSeg);
    for kt = 1:length(cID)
        cnt = cnt+1;
        selInfo(cnt,:) = [it cID(kt)];
        ts = spk.SpikeTimesSeg( spk.assignedClusterSeg==cID(kt)).*1e3;
        dx = diff(ts);
        dx(sign(dx)==-1) = [];
        n = hist(dx,dt);
        isiH{cnt} = n(1:end-1);
        wvf{cnt} = spk.wavf( spk.oriIx( spk.assignedClusterSeg==cID(kt) ),:);
        wC{cnt} = wltCoeffs{it}(spk.oriIx(spk.assignedClusterSeg==cID(kt)),:);
    end;
    
end;

%%
C = [1 0 0;0 1 0;0 0 1; 0 0 0];
mark = {'x','o','^','d','s'};

selIx = [];
mwID = unique(selInfo(:,1));
figure;
for it = 1:length(mwID)
    
    ix = find(selInfo(:,1)==mwID(it));
    for jt = 1:length(ix)
        
        k = length(ix);
        if k==1
            k = 2;
        end;
        
        subplot(4,k,jt);
        hold on;
        plot(linspace(0,2,64),wvf{ix(jt)}(1:3:end,:),'Color',C(jt,:));
        plot(linspace(0,2,64),mean(wvf{ix(jt)},1),'k');
        axis tight;
        
        subplot(4,k,jt+length(ix));
        bar(0:250,isiH{ix(jt)});
        xlim([-10 250]);                
        
        subplot(4,k,[3*k-1 4*k-1]);
        hold on;
        plot(wC{ix(jt)}(:,1),wC{ix(jt)}(:,2),mark{jt},'Color',C(jt,:),'MarkerSize',8);
        axis tight;
        
        subplot(4,k,[3*k 4*k]);
        hold on;
        plot(linspace(0,2,64),wvf{ix(jt)}(1:3:end,:),'Color',C(jt,:));
        plot(linspace(0,2,64),mean(wvf{ix(jt)},1),'k');        
    end;
    set(gcf,'Color','w');
    
    for jt = 1:length(ix)
        [s] = input(['unit[1],noise[2],',num2str(jt),'/',num2str(length(ix))],'s');
        if strcmp(s,'1');
            selIx = [selIx ix(jt)];
        end;
    end;
    clf;
    
end;

%%
x = zeros(1,20);
x(15) = 1;

figure;
plot(x);

Tau = 1:20;
t = 1;
expKern = 1-exp(-t./Tau);

y = filter(expKern,1,spkTrains(3,:));
figure;
plot(y);

%%
dt = -500:5000;
eR = [];
for it = 1:length(trlENC)
    
    spkTrains = [];
    ensResp = {};
    bflab = {};
    cnt = 0;
    for jt = 1:length(sortedSpikes)
    
        [spk] = sortedSpikes{jt};
        
        cID = unique(spk.assignedClusterSeg);
        for kt = 1:length(cID)            
            [ts] = spk.SpikeTimesSeg(spk.trl == trlENC(it) & spk.assignedClusterSeg == cID(kt)).*1e3;
            [ts] = ts(ts>=-5e2 & ts<=5e3);
            %if ~isempty(ts)
                cnt = cnt+1;
                ensResp{cnt} = ts;
                bflab{cnt} = BFlab{jt};
                spkTrains(cnt,:) = hist(ts,dt);
            %;    
        end;
        
    end;
    ensResp = ensResp(selIx);
    bflab = bflab(selIx);
    hemLab = cell(1,length(bflab));
    
    for jt = 1:length(bflab)
        hemLab(jt) = {bflab{jt}(end)};
    end;
    [~,sIx] = sort(hemLab);
    bflab = bflab(sIx);
    ensResp = ensResp(sIx);
    
    for jt = 1:length( ensResp )
        eR(it,jt,:) = hist( ensResp{jt},dt );
    end;
    
%     figure;
%     hold on;
%     plot([0 0],[0 length(ensResp)+1],'Color',[.75 .75 .75]);
%     plot([2e3 2e3],[0 length(ensResp)+1],'Color',[.75 .75 .75]);
%     for jt = 1:length( ensResp )
%         
%         st = ensResp{jt};
%         st = st(st>=-500 & st <=5e3);
%         
%         x = [st;st];
%         y = jt*ones(1,length(x));
%         y = [y-.5;y+.5];
%         line(x,y,'Color',cIx{strcmp(BFid,bflab(jt))});
%     end;
%     xlim([-500 5000]);
%     set(gca,'YTick',1:2:length(ensResp));
%     set(gca,'YTickLabel',bflab(1:2:length(ensResp)));
    
end;

%%
sm = [];
for it = 1:size(eR,1)
    for jt = 1:size(eR,2)
        
        x  = squeeze(eR(it,jt,:));
        
        %%
        T = 20;
        tf = find(x==1);
        e = zeros(1,length(x));
        for kt = 1:length(tf)
            for t = 1:length(x)
                if t<tf(kt)
                else
                    e(t) = e(t)+exp(-(t-tf(kt))/T);
                end;
            end;
        end;
        %%
        sm(it,jt,:) = e;
    end;
end;

%%
[cueDat] = sm(:,:,dt>=0 & dt<=2e3);
[encDat] = sm(:,:,dt>=2e3 & dt<=4e3);

X = [cueDat;encDat];
Y = [ones(1,size(cueDat,1)) 2*ones(1,size(encDat,1))]';

rIx = randperm(size(X,1));
trIx = rIx(1:fix(length(rIx)*0.8));
rIx = setdiff(rIx,trIx);
rIx = rIx(randperm(length(rIx)));
cvIx = rIx(1:fix(length(rIx)*0.5));
teIx = setdiff(rIx,cvIx);
rIx = setdiff(rIx,[teIx cvIx trIx]);

err = [];
for it = 1:size(cueDat,3)
        
    OBJ=fitcecoc(squeeze(X(trIx,:,it)),Y(trIx));
    Yp = predict(OBJ,squeeze(X(cvIx,:,it)));
    
    err(it) = sum(abs(diff([Yp Y(cvIx)],[],2)))/length(cvIx);
    
end;
