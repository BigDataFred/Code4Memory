function [selInfo,selIx] = VisualizeAndScreeningOfClusterData(spkDat)

cIx = {[255 0 0]/255,[51 255 51]/255,[0 0 204]/255,[255 128 0]/255,[255 255 0]/255,[0 255 255]/255};

%%
BFlab = cell(1,length(spkDat.chanLab));
for it = 1:length( spkDat.chanLab)
    
    x = spkDat.chanLab{it};
    x(end) = [];
    BFlab(it) ={x};
end;
BFid = unique(BFlab)';

%%
cnt = 0;
wvf = {};
isi = {};
dt = 0:251;
selInfo = [];
wC = {};
spkCnt = [];
for it = 1:length(spkDat.sortedSpikes)
    it
    [spk] = spkDat.sortedSpikes{it};
    
    cID = unique(spk.assignedClusterSeg);
    for kt = 1:length(cID)
        cnt = cnt+1;
        selInfo(cnt,:) = [it cID(kt)];
        ts = spk.SpikeTimesSeg( spk.assignedClusterSeg==cID(kt)).*1e3;
        spkCnt(cnt) = length(ts);
        dx = diff(ts);
        dx(sign(dx)==-1) = [];
        n = hist(dx,dt);
        isiH{cnt} = n(1:end-1);
        wvf{cnt} = spk.wavf( spk.oriIx( spk.assignedClusterSeg==cID(kt) ),:);
        wC{cnt} = spkDat.wltCoeffs{it}(spk.oriIx(spk.assignedClusterSeg==cID(kt)),:);
    end;
    
end;

%%
C = [1 0 0;0 1 0;0 0 1; 0 0 0; 1 1 0; 0 1 1;.75 .75 .75];
mark = {'x','o','^','d','s'};

selIx = [];
mem = [];

mwID = unique(selInfo(:,1));
figure;
it = 1;
while it <= length(mwID)
    
    f = 0;
    while f<1
        ix = find(selInfo(:,1)==mwID(it));
        for jt = 1:length(ix)
            
            k = length(ix);
            if k==1
                k = 2;
            end;
            
            subplot(4,k,jt);
            hold on;
            plot(linspace(0,2,64),wvf{ix(jt)}(1:5:end,:),'Color',C(jt,:));
            plot(linspace(0,2,64),mean(wvf{ix(jt)},1),'k');
            axis tight;
            title(spkCnt(ix(jt)));
            
            subplot(4,k,jt+length(ix));
            bar(0:250,isiH{ix(jt)});
            xlim([-10 250]);
            
            subplot(4,k,[3*k-1 4*k-1]);
            hold on;
            plot(wC{ix(jt)}(:,1),wC{ix(jt)}(:,2),mark{jt},'Color',C(jt,:),'MarkerSize',8);
            axis tight;
            
            subplot(4,k,[3*k 4*k]);
            hold on;
            plot(linspace(0,2,64),wvf{ix(jt)}(1:5:end,:),'Color',C(jt,:));
            plot(linspace(0,2,64),mean(wvf{ix(jt)},1),'k');
            
        end;
        set(gcf,'Color','w');
        
        a = [];
        for jt = 1:length(ix)
            [s] = input(['unit[1],noise[2],',num2str(jt),'/',num2str(length(ix))],'s');
            if strcmp(s,'1');
                selIx = [selIx ix(jt)];
                a = [a ix(jt)];               
                f = 1;
            elseif strcmp(s,'r');
                it = it-1;
                if ~isempty(mem)
                    selIx(end-(mem(end)-1):end) = [];
                    selInfo(selIx(end-(mem(end)-1):end),3) = 0;
                end;
                clf;
                break;
            else
                f = 1;
            end;
        end;
        
        if length(a) > 1
            mem = [mem length(a)];
            [s] = input(['enter idx to merge clusters [1,2,3;..,n]'],'s');
            s = str2num(s);
            if ~isempty(s)
                for kt = 1:size(s,1)
                    for lt = 1:size(s,2)
                        selInfo(ix(s(kt,lt)),3) = kt;
                    end;
                end;
            end;
        end;
    end;   
    it = it+1;
    clf;
    
end;