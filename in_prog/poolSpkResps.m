%%
addpath('/home/rouxf/prj/Bham/code/mcode/utils');

%%
pID = {'P02','P04','P05','P07','P08'};

FR1 = [];
FR2 = [];
frH1 = [];
frM1 = [];
frH2 = [];
frM2 = [];
acfB = [];
acfC = [];
acfE = [];
dt  = [-1e3:5e3];
dt2 = [-1e3:250:5e3];

%%
for pIt=1:length(pID)
	
	p2dOrig=['/media/rouxf/rds-share/iEEG_DATA/MICRO/',pID{pIt},'/spkRes/'];
	spkF = dir([p2dOrig,pID{pIt},'_*_*_spkResp.mat']);
	
	for it = 1:length(spkF)
		fprintf([spkF(it).name,'\n']);
		spkDat = load([p2dOrig,spkF(it).name])
		
		[hitIdx] = find(ismember(spkDat.hitIdx,sort([spkDat.hitIdx;spkDat.missIdx])));
		[missIdx] = find(ismember(spkDat.missIdx,sort([spkDat.hitIdx;spkDat.missIdx])));

        fr1 = [];
		fr2 = [];	
		spkTs = [];
		for kt = 1:length(spkDat.sigSel)
			isiH = spkDat.isiH{spkDat.sigSel(kt)};
			pct = (sum(isiH(1:4))/sum(isiH))*100;
			if pct <3
				
				[spk] = spkDat.spkRaster{spkDat.sigSel(kt)};
				x = [];
				x(1,1) = [mean(sum(spk(hitIdx,dt>=2e3&dt<4e3),2)/2)];
				x(1,2) = [mean(sum(spk(hitIdx,dt>=0&dt<2e3),2)/2)];
				x(1,3) = [mean(sum(spk(missIdx,dt>=2e3&dt<4e3),2)/2)];                        	
                x(1,4) = [mean(sum(spk(missIdx,dt>=0&dt<2e3),2)/2)];				
                
				x2 = spkDat.n2{spkDat.sigSel(kt)};
				x2 = sum(x2,1)/size(x2,1)/0.25;
				x2 = (x2-mean(x2(dt2>=-1e3&dt2<=0)))/std(x2(dt2>=-1e3&dt2<=0));
				[~,mIx] = max(x2); 
                if dt2(mIx)<2e3
                    FR1 = [FR1;x2];
                    fr1 = [fr1;x(1,[2 4])];
                else
                    FR2 = [FR2;x2];
                    fr2 = [fr2;x(1,[1 3])];
                end;
                
				ix = [];for it = 1:size(spk,1);ix(it) = size(spk,1)*(it-1)+(it-1)+1;end;	
				[xc1] = xcorr(spk(:,dt>=0 & dt<=2000)',500);
				[xc2] = xcorr(spk(:,dt>=2000 & dt<=4000)',500);
                [xc3] = xcorr(spk(:,dt>=-2000 & dt<=0)',500);
                acfB = [acfC;mean(xc3(:,ix),2)'];
				acfC = [acfC;mean(xc1(:,ix),2)'];
				acfE = [acfE;mean(xc2(:,ix),2)'];	

				spkTs = [spkTs;spk];	
			end;
            %fr2 = zscore(fr2(:));
		%fr2 = reshape(fr2,[length(fr2)/4 4]);
        if ~isempty(fr1)
            frH1 = [frH1;fr1(:,1)];
            frM1 = [frM1; fr1(:,2)];
        end;
        
        if ~isempty(fr2)
            frH2 = [frH2;fr2(:,1)];
            frM2 = [frM2; fr2(:,2)];
        end;
        
		end;
		clear fr2 fr1;
	end;
end;
fprintf(['Number of units included :',num2str(size(FR2,1)),'\n']);

%%
figure;
a = [];
for it = 92%1:size(acfC,1);
    subplot(131);
    a = [a gca];
    x = acfB(it,:);
    x(501) = NaN;
    bar(-500:500,x,'FaceColor',[.75 .75 .75]);
    subplot(132);
    a = [a gca];
    x = acfC(it,:);
    x(501) = NaN;
    bar(-500:500,x,'FaceColor',[.75 .75 .75]);
    subplot(133);
    a = [a gca];
    x = acfE(it,:);
    x(501) = NaN;
    bar(-500:500,x,'FaceColor',[.75 .75 .75]);
end;
axis tight;
for it = 1:length(a)
    ylabel(a(it),'Coincidences');
    xlabel(a(it),'Lag (ms)');
end;
set(a,'YTick',[0 3]);
set(a,'XTick',[-250:125:250]);
set(a,'XTickLabel',{'-250','','','','250'});
set(a,'XLim',[-250 250]);
set(a,'YLim',[0 3]);
set(a,'LineWidth',2);
set(a,'FontName','Arial');
set(gcf,'Color','w');

%%
M1 = mean(FR1,1);
M1 = conv(M1,gausswin(5),'same')./sum(gausswin(5));
SE1 = conv(std(FR1,0,1)./sqrt(size(FR1,1)),gausswin(5),'same')./sum(gausswin(5));
M2 = conv(mean(FR2,1),gausswin(5),'same')./sum(gausswin(5));
SE2 = conv(std(FR2,0,1)./sqrt(size(FR2,1)),gausswin(5),'same')./sum(gausswin(5));

figure;
hold on;
jbfill(dt2,M1-SE1,M1+SE1,[0 128 0]./255,[0 128 0]./255,1,.4);
hold on;
jbfill(dt2,M2-SE2,M2+SE2,[255 140 0]./255,[255 140 0]./255,1,.4);
hold on;
plot(dt2,M1,'Color',[0 128 0]./255,'LineWidth',2,'LineSmoothing','on');
plot(dt2,M2,'Color',[255 140 0]./255,'LineWidth',2,'LineSmoothing','on');

plot([min(dt2) max(dt2)],[0 0],'k--','LineWidth',2,'LineSmoothing','on');
set(gca,'LineWidth',2);
set(gca,'XTick',[0:2e3:4e3]);
set(gca,'YTick',[0 2]);
xlabel('Time rel. to cue onset (ms)','Fontsize',14);
ylabel('Norm. firing rate (\sigma)','Fontsize',14);
xlim([-500 4750]);

figure;
hold on;
m1 = mean(FR1(:,dt2 >=0 & dt2 <=2e3),2);
m2 = mean(FR1(:,dt2 >=2e3 & dt2 <=4e3),2);
m3 = mean(FR2(:,dt2 >=0 & dt2 <=2e3),2);
m4 = mean(FR2(:,dt2 >=2e3 & dt2 <=4e3),2);
sd1 = std(m1)/sqrt(length(m1));
sd2 = std(m2)/sqrt(length(m1));
sd3 = std(m3)/sqrt(length(m1));
sd4 = std(m4)/sqrt(length(m1));

errorbar(1:2,mean([m1 m2],1),[sd1 sd2],'s-');
errorbar(1:2,mean([m3 m4],1),[sd3 sd4],'o-');

%%
figure;
x = [frH1(:,1) frM1(:,1)];
[delIx,~] = ind2sub(size(x),find(isnan(x)));
x(delIx,:) = [];
frH1(delIx) = [];
frM1(delIx) = [];
M = mean(x,1);
SE = std(x,0,1)/sqrt(size(x,1)-1);

subplot(121);
a = gca;
hold on;
for it = 1:size(x,1)
	plot([1 2],[frH1(it,1) frM1(it,1)],'-','Color',[.75 .75 .75]);
end;
plot(ones(size(x,1),1),frH1(:,1),'bo');
plot(2*ones(size(x,1),1),frM1(:,1),'ro');
title('Cue');

x = [frH2(:,1) frM2(:,1)];
[delIx,~] = ind2sub(size(x),find(isnan(x)));
x(delIx,:) = [];
frH2(delIx) = [];
frM2(delIx) = [];
M = mean(x,1);
SE = std(x,0,1)/sqrt(size(x,1)-1);
subplot(122);
a = [a gca];
hold on;
for it = 1:size(x,1)
	plot([1 2],[frH2(it,1) frM2(it,1)],'-','Color',[.75 .75 .75]);
end;
plot(ones(size(x,1),1),frH2(:,1),'bo');
plot(2*ones(size(x,1),1),frM2(:,1),'ro');
title('Encoding');

set(a,'XLim',[0 3]);
set(a,'XTick',1:2);
set(a,'XTickLabel',{'LR','LF'});
set(a,'LineWidth',2);
for it = 1:length(a)
	ylabel(a(it),'Firing rate (Spikes/s)');
end;

