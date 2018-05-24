%%
addpath('/home/rouxf/prj/Bham/code/mcode/params/');
addpath(genpath('/home/rouxf/tbx/osort-v3-rel/code/'));
addpath(genpath('/home/rouxf/tbx/chronux_2_11/'));

addpath(genpath('/home/rouxf/prj/Bham/code/mcode/utils/'));

Fs =  32e3;
[par] = set_parameters_Bham(Fs);

%%
pID = 'P07';

rpath = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',pID,'/'];

chck = dir(rpath);
chck = chck([chck.isdir] == 1);
ix = regexp({chck.name},'EM');
sel=[];c=0;for it = 1:length(ix); if ~isempty(ix{it});c = c+1;sel(c) = it;end;end;
expMode = {chck(sel).name};

%%
for expSel = 1:length( expMode )
    
    sesh = dir([rpath,expMode{expSel},'/']);
    chck = regexp({sesh(:).name},'\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}')';
    sel=[];c=0;for it = 1:length(chck); if ~isempty(chck{it});c = c+1;sel(c) = it;end;end;
    sesh = sesh(sel);
    sesh = sesh([sesh(:).isdir]==1);
    sesh = {sesh.name};
    
    for seshSel = 1%:length( sesh )
        
        path2file = [rpath,expMode{expSel},'/',sesh{seshSel},'/'];
        savepath = path2file;
        fn = [pID,'_',expMode{expSel},'_',sesh{seshSel},'_spkDataStimLockedSegmented.mat'];        
        load([path2file,fn]);
        
        %% get the trial events
        trl = [];
        for it = 1:length( sortedSpikes )
            trl = [trl sortedSpikes{it}.trl];
        end;
        trl = unique( trl );
        
        %% sort the spike time data from left to right
        chck = regexp(chanLab,'L');
        sel = [];
        for it = 1:length(chck)
            if ~isempty(chck{it})
                sel = [sel it];
            end;
        end;
        ixL = sel;
        ixR = setdiff(1:length(chanLab),ixL);
        ix = [ixL ixR];
        sortedSpikes = sortedSpikes(ix);
        chanLab = chanLab(ix);
        wltCoeffs = wltCoeffs(ix);
        
        %% visualize the spike time data
        flag = 0;
        %figure;
        while flag <1
            
            uSel = {};
            for ut = 1:3
                c = 0;
                for it = 1:length( sortedSpikes )
                    
                    cID = unique( sortedSpikes{it}.assignedClusterSeg);
                    
                    for jt = 1:length( cID )
                        
                        c = c+1;
                        
                        [ix] = find( sortedSpikes{it}.assignedClusterSeg == cID( jt ) );
                        
                         wvf  = sortedSpikes{it}.wavf(sortedSpikes{it}.oriIx( ix ),:);
                        
                        [~,mx] = max(abs(mean(wvf,1)));
                        ref = sign(mean(wvf(:,mx),1));
                        
                        sn = [];
                        for zt = 1:size(wvf,1);
                            sn(zt) = sign(wvf(zt,mx));
                        end
                        delIx = find(sn ~=ref);
                        wvf(delIx,:) = [];
                        ix(delIx) = [];
                        
                        ts = sortedSpikes{it}.SpikeTimesSeg(ix).*1e3;
                        trl2 = sortedSpikes{it}.trl(ix);
                        
                        dx = diff(ts);
                        dx(dx<0) = [];
                        
                        pct = length(find(dx < 3))/length(dx)*100;
                        dx(dx>250) = [];
                        [isiH,~] = hist(dx,[0:250]);
                        
                        FR = [];
                        n = [];
                        xc = [];
                        n2 = [];
                        x = {};
                        dt = -1000:1:5e3;
                        dt2 = -500:250:5e3;
                        for kt = 1:length(trl)
                            
                            x{kt} = ts(trl2 == trl(kt));
                            x{kt} = x{kt}(x{kt} > dt2(1) & x{kt} <= dt2(end));
                            
                            [n(kt,:),~] = hist(x{kt},dt2);
                            
                            FR(kt) = length(x{kt})/((dt(end)-dt(1))./1e3);
                            
                            [n2(kt,:),~] = hist(x{kt},dt);                            
                            [xc(kt,:),lag] = xcorr(n2(kt,:),250);
                            
                        end;
                        
%                         lag = 300;
%                         XC2 = zeros(length(dt)-(lag/2+1),2*lag+1);
%                         
%                         for mt = 1:size(n2,1)
%                             xc = [];
%                             parfor kt = 1:length(dt)
%                                 if (kt > lag/2) && (kt < (length(dt)-lag/2))
%                                     xc(kt,:) = xcorr(n2(mt,kt-lag/2:kt+lag/2),lag);
%                                 end;
%                             end;
%                             XC2 = XC2+xc;
%                         end;
%                         XC2 = XC2./mt;
%                         XC2(:,lag+1) = NaN;
%                         
%                         params = [];
%                         params.Fs = 1e3;
%                         params.pad = 0;
%                         params.tapers = [25 49];
%                         
%                         [S,f,R] = mtspectrumpb(n2',params);                                                                      
                        
                        wvfd = [];
                        %dx2 = linspace(min(min(wvf)),max(max(wvf)),200);
                        dx2 = linspace(min(mean(wvf,1)-2*max(std(wvf,0,1))),max(mean(wvf,1)+2*max(std(wvf,0,1))),200);
                        for kt = 1:size(wvf,2)
                            x2  = wvf(:,kt);
                            [nx,~] = hist(x2,dx2);
                            nx = nx./sum(nx);
                            wvfd(kt,:) = nx;
                        end;
                        
                        SNR = max(abs(mean(sortedSpikes{it}.wavf(sortedSpikes{it}.oriIx( ix ),:),1)))/sortedSpikes{it}.SD;
                        SNR = round(SNR*10)/10;
                        
                        if any(dx<0)
                            error('toto');
                        end;
                        
                        figure
                        subplot(3,2,1);
                        hold on;
                        imagesc( linspace(0,2,64), dx2 , wvfd' );
                        plot(linspace(0,2,64),mean(mean(sortedSpikes{it}.wavf(sortedSpikes{it}.oriIx( ix ),1:10)))*ones(1,64),'w')
                        axis xy;axis tight;%caxis([0 .1]);
                        title([chanLab{it},' cID:', num2str(jt)]);
                        xlabel('Time [ms]');
                        ylabel('Amolitude [\muV]');
                        
                        subplot(3,2,2);
                        hold on;
                        plot( linspace(0,2,64), wvf ,'r');
                        h = plot( linspace(0,2,64), mean( wvf ,1) ,'k');
                        axis tight;
                        legend(h,['SNR: ',num2str(SNR)]);
                        legend boxoff;
                        title(['Number of spikes:', num2str(length(ix))]);
                        xlabel('Time [ms]');
                        ylabel('Amplitude [\muV]');
                        
                        subplot(3,2,3);
                        bar(0:250,isiH);
                        axis tight;xlim([-10 250]);
                        title([num2str(round(pct*100)/100),'pct < 3ms']);
                        xlabel('Time [ms]');
                        ylabel('Count');
                        
                        subplot(3,2,4);
                        M = mean(xc(:,(length(lag)-1)/2+2:end),1);
                        SE = std(xc(:,(length(lag)-1)/2+2:end),0,1)./sqrt(size(xc,1)-1);
                        [ph,msg] = jbfill(lag((length(lag)-1)/2+2:end),(M-SE),(M+SE),[0 1 1],[0 1 1],0,1);
                        hold on;
                        plot(lag((length(lag)-1)/2+2:end),M,'b');
                        axis tight;
                        xlabel('Lag [ms]');
                        ylabel('Coincidences');
                        title('Mean Autocorrelation');
                        if max(get(gca,'YLim')) < 1
                            ylim([0 1]);
                        end;
                        xlim([-10 max(lag)]);
                        
                        subplot(3,2,5);
                        hold on;
                        
                        for kt = 1:length(trl)
                            y = kt*ones(1,length(x{kt}));
                            y = [y-.5;y+.5];
                            line([x{kt};x{kt}],y,'Color','k');
                        end;
                        plot([0 0],[0 length(trl)+1],'r');
                        plot([2e3 2e3],[0 length(trl)+1],'r');
                        xlim([-500 5e3]);ylim([0 length(trl)+1]);
                        title('Rastergram');
                        xlabel('Time rel. to Cue Onset [ms]');
                        ylabel('Trial event #');
                        
                        subplot(3,2,6);
                        hold on;
                        M = sum(n,1)/size(n,1)/2.5e-1;
                        SE = 1/size(n,1)*sum((n-ones(size(n,1),1)*M).^2,1);
                        SE = SE/sqrt(size(n,1)-1);
                        plot(dt2,M,'ks-','MarkerFaceColor','k');
                        for kt = 1:length(dt2)
                            plot(dt2(kt)*ones(1,2),[M(kt) M(kt)+SE(kt)],'k');
                            plot([dt2(kt)-50 dt2(kt)+50],M(kt)+SE(kt)*ones(1,2),'k');
                            plot([dt2(kt)-50 dt2(kt)+50],M(kt)+SE(kt)*ones(1,2),'k');
                        end;
                        plot([0 0],[0 max(M+SE)],'r');
                        plot([2e3 2e3],[0 max(M+SE)],'r');
                        axis tight;xlim([-500 5e3]);
                        xlabel('Time rel. to Cue Onset [ms]');
                        ylabel('Spike/s');
                        title(['Average Firing Rate:', num2str(round(mean(FR)*10)/10),'Hz']);
                        
%                         [s] = input('enter SUA [1], MUA [2] or noise [3]','s');
%                         uSel{ut}(c,:) = [it jt str2double(s)];
                        
                    end;
                    
                end;
            end;
            
%             x = [];
%             for it = 1:length( uSel )
%                 x = [x uSel{it}(:,3)];
%             end;
%             rxy = corr(x);
%             rxy = [rxy(2:3,1)' rxy(3,2)];
%             rxy
%             
%             if ( sum(rxy >= 0.97)==3 )
%                 flag = 1;
%             end;
            
        end;
    end  
end
%         %%        
%         save([savepath,pID,'_',expMode{expSel},'_',sesh{seshSel},'_unitClassification.mat'],'uSel','rxy');
%         
%         %%
%         [mRxy,mIx] = max(rxy)
% 
%         if mIx ==1% 1->2
%             uSel = uSel{2};
%         elseif mIx ==2% 1->3
%             uSel = uSel{3};
%         else% 2->3
%             uSel = uSel{3};
%         end;
%         
%         %% remove noise
%         delIx = find(uSel(:,3) ==3);
%         uSel(delIx,:) = [];
%         
%         %%
%         BFlab = uSel(:,1);
%         BFid  =unique(uSel(:,1));
%         
%         uSel2 = zeros(size(uSel));
%         uSel2 = [uSel2 zeros(size(uSel,1),1)];
%         for it = 1:length(BFid);
%             ix = find(uSel(:,1) == BFid(it));
%             x = uSel(ix ,:);
%             u = x(:,3)==2;
%             if sum(u)>1
%                 x = [x u];
%             else
%                 x = [x zeros(length(u),1)];
%             end;
%             uSel2(ix,:) = x;
%         end;
% 
%         %%
%         %BFlab = cell(size(uSel,1),1);
%         %for it = 1:size(uSel,1)
%         %    
%         %    BFlab(it,1) = {chanLab{uSel(it,1)}(1:regexp(chanLab{uSel(it,1)},'\d{1}')-1)};
%         %    
%         %end;
%         %BFid = unique(BFlab);
%         
%         BFlab = uSel(:,1);
%         BFid  =unique(uSel(:,1));
%         
%         for it = 1:length( BFid )
%             
%             %selIx = find(strcmp(BFlab,BFid(it)) & uSel(:,3)==1);
%             selIx = find((BFlab == BFid(it)) & uSel(:,3)==1);
%             
%             lab = uSel(selIx,1);
%             
%             if length(lab) >1
%                 
%                 c = 0;
%                 p = [];
%                 d=[];
%                 for kt = 1:length(lab)
%                     for lt = kt+1:length(lab)
%                         c = c+1;
%                         p(c,:) = [kt lt];
%                     end;
%                 end;
%                 
%                 c=0;
%                 for kt = 1:size(p,1)
%                     c = c+1;
%                     d1 = uSel(selIx(p(kt,1)),:);
%                     d2 = uSel(selIx(p(kt,2)),:);
%                     
%                     ix1 = find( sortedSpikes{d1(1)}.assignedClusterSeg == d1(2) );
%                     ix2 = find( sortedSpikes{d2(1)}.assignedClusterSeg == d2(2) );
%                     
%                     ts1 = sortedSpikes{d1(1)}.SpikeTimesSeg(ix1).*1e3;
%                     ts2 = sortedSpikes{d2(1)}.SpikeTimesSeg(ix2).*1e3;
%                     ts1= ts1(ts1 >= -500 & ts1 <= 5e3);
%                     ts2= ts2(ts2 >= -500 & ts2 <= 5e3);
%                     
%                     dt = -500:5e3;
%                     [n1,~] = hist(ts1,dt);
%                     [n2,~] = hist(ts2,dt);
%                     [xc,lag] = xcorr(n1,n2,250);
%                                                             
%                     spikes1 = sortedSpikes{d1(1)}.wavf_decorr(sortedSpikes{d1(1)}.oriIx(ix1),:);
%                     spikes2 = sortedSpikes{d2(1)}.wavf_decorr(sortedSpikes{d2(1)}.oriIx(ix2),:);                   
%                     
%                     C1 = wltCoeffs{d1(1)}(sortedSpikes{d1(1)}.oriIx(ix1),:);
%                     C2 = wltCoeffs{d2(1)}(sortedSpikes{d2(1)}.oriIx(ix2),:);
%                     
%                     cf1 = randperm(size(C1,1));
%                     teIx1 = cf1(1:floor(size(C1,1)*0.1));%10%
%                     cf1(1:length(teIx1)) = [];
%                     teIx3 = cf1(1:floor(size(C1,1)*0.1));%10%
%                     
%                     cf2 = randperm(size(C2,1));
%                     teIx2 = cf2(1:floor(size(C2,1)*0.1));% 10%                
%                     cf2(1:length(teIx2)) = [];
%                     teIx4 = cf2(1:floor(size(C2,1)*0.1));% 10%
%                     
%                     trIx1 = 1:size(C1,1);
%                     trIx2 = 1:size(C2,1);
%                     trIx1([teIx1 teIx3]) = [];% 80 %
%                     trIx2([teIx2 teIx4]) = [];% 80 %
%                     
%                     p2 = [];
%                     c2 = 0;
%                     for lt = 1:size(C1,2)
%                         for mt = lt+1:size(C1,2)
%                             c2 = c2+1;
%                             p2(c2,:) = [lt mt];
%                         end;
%                     end;
%                     
%                     acc = zeros(size(p2,1),3);
%                     for lt = 1:size(p2,1)                                                
%                         X = [C1(trIx1,p2(lt,:));C2(trIx2,p2(lt,:))];
%                         Y = [zeros(length(trIx1),1);ones(length(trIx2),1)];                        
%                         LDAdsc = fitcdiscr(X,Y);
%                         
%                         cM = confusionmat(Y,predict(LDAdsc,X));                        
%                         accTR = sum([cM(1,1) cM(2,2)])/sum(cM(:));
%                         
%                         X = [C1(teIx1,p2(lt,:));C2(teIx2,p2(lt,:))];
%                         Y = [zeros(length(teIx1),1);ones(length(teIx2),1)];                        
%                         cM = confusionmat(Y,predict(LDAdsc,X)); 
%                         accCV = sum([cM(1,1) cM(2,2)])/sum(cM(:));
%                         
%                         X = [C1(teIx3,p2(lt,:));C2(teIx4,p2(lt,:))];
%                         Y = [zeros(length(teIx3),1);ones(length(teIx4),1)];                        
%                         cM = confusionmat(Y,predict(LDAdsc,X)); 
%                         accTE = sum([cM(1,1) cM(2,2)])/sum(cM(:));
%                         
%                         acc(lt,:) = [accTR accCV accTE]; % training, cross-val, test                       
%                     end;
%                     
%                     [~,ix] = max(acc(:,3));
%                     acc = acc(ix,:);
%                     cix = p2(ix,:);
%                     
%                     [m1,m2, projectedResidual1,projectedResidual2,overlap,d(c)] = projectionTest( spikes1,spikes2 );
%                     [bc1,fhat1,h1]=estimatePDF(projectedResidual1,30);
%                     [bc2,fhat2,h2]=estimatePDF(projectedResidual2,30);
%                     
%                     %this is used for calculating Rsquare
%                     dist1 = normpdf(bc1,0,1);
%                     dist2 = normpdf(bc2-d(c),0,1);
%                     [Rsquare1] = calcRsquare( fhat1, dist1 );
%                     [Rsquare2] = calcRsquare( fhat2, dist2 );
%                     
%                     %more sampling points for plotting purposes
%                     xDist=-4:.1:4;
%                     distPlot=normpdf(xDist,0,1);
%                     
%                     figure;
%                     subplot(3,2,1);
%                     hold on;                    
%                     ha1 = plot(linspace(0,2,size(spikes1,2)),spikes1,'r');
%                     ha2 = plot(linspace(0,2,size(spikes2,2)),spikes2,'b');
%                     axis tight;
%                     xlabel('Time [ms]');
%                     ylabel('Amplitude [\muV]');
%                     title(chanLab{d1(1)});
%                     legend([ha1(1);ha2(1)],['clusterID:',num2str(d1(2))],['clusterID:',num2str(d2(2))]);
%                     subplot(3,2,3);
%                     hold on;                    
%                     plot(linspace(0,2,size(spikes1,2)),mean( spikes1 ,1),'r','LineWidth',3);
%                     plot(linspace(0,2,size(spikes2,2)),mean( spikes2 ,1),'b','LineWidth',3);                    
%                     axis tight;
%                     xlabel('Time [ms]');
%                     ylabel('Amplitude [\muV]');
%                     title(chanLab{d2(1)});
%                     legend(['clusterID:',num2str(d1(2))],['clusterID:',num2str(d2(2))]);
%                     
%                     subplot(3,2,2);
%                     hold on;                    
%                     plot(lag,xc,'k');
%                     [~,mIx] = max(xc);
%                     plot(lag(mIx),xc(mIx),'ro');                   
%                     axis tight;
%                     
%                     title(['Max: ',num2str(lag(mIx)),'ms']);
%                     xlabel('Lag [ms]');
%                     ylabel('Coincidences');
%                     
%                     subplot(3,2,[4 6]);
%                     bar(bc1, fhat1, 'r');
%                     hold on
%                     h=bar(bc2, fhat2, 'b');
%                     h = [];
%                     h(1) = plot(xDist,distPlot,'r','LineWidth',2.5);
%                     h(2) = plot(xDist+d(c),distPlot,'b','LineWidth',2.5);
%                     hold off;
%                     axis tight;
%                     ylim([0 0.5]);
%                     title(['d =',num2str(d(c))]);
%                     legend(h,['R^2:',num2str(round(Rsquare1*100)/100)],['R^2:',num2str(round(Rsquare2*100)/100)],'location','NorthWest');
%                     legend boxoff;
%                     
%                     subplot(3,2,5);
%                     hold on;
%                     X = wltCoeffs{d1(1)}(:,cix(1));
%                     Y = wltCoeffs{d1(1)}(:,cix(2));
%                     plot(X,Y,'.','Color',[.75 .75 .75]);
%                     plot(C1(:,cix(1)),C1(:,cix(2)),'o','Color',[.9 0 0]);
%                     plot(C2(:,cix(1)),C2(:,cix(2)),'x','Color',[0 0 .9]);
%                     axis tight;                    
%                     sx = [];
%                     for ht = 1:length(acc)
%                         sx = [sx,' ',num2str(round(acc(ht)*1e2)/1e2)];
%                     end;
%                     title(['Classification Accuracy (%):',sx]);
%                     xlabel('X1');
%                     ylabel('X2');
%                     
%                     [s] = input('Do you want to merge these clusters? [Y/N]','s');
%                     if strcmp(s,'Y') || strcmp(s,'y')
%                         uSel2(selIx,end) = ones(length(selIx),1);
%                     end;
%                 end;
%                 
%             end;
%         end;
% 
%     end;
% end;
% 
% %% 
% [BFid] = unique(uSel2(:,1));
% dt = -1500:250:7000;
% dt2 = -1500:7000;
% 
% for it = 1:length( BFid )
%     
%     
%     ix = find( (uSel2(:,1) == BFid(it)) & (uSel2(:,3) ==2) );
%     
%     if length(ix) >1
%         m = find( uSel2(ix,end) );
%     else
%         m = uSel2(ix,2);
%     end;
%     
%     cIx = ismember(sortedSpikes{BFid(it)}.assignedClusterSeg,m);
%     
%     if any(cIx ~=0)
%         trl = sortedSpikes{BFid(it)}.trl(cIx);
%         ts = sortedSpikes{BFid(it)}.SpikeTimesSeg(cIx).*1e3;
%         
%         trl(ts<-1500) = [];
%         ts(ts<-1500) = [];
%         
%         dx = diff(ts);
%         dx(sign(dx)==-1) = [];
%         [isi] = hist(dx,0:201);
%         isi(end) = [];
%         
%         trlID = unique(trl);
%         
%         figure;
%         subplot(6,4,[1 2 5 6 9 10 13 14 17 18]);
%         hold on;
%         n = [];
%         xc = [];
%         for jt = 1:length( trlID )
%             
%             trlIx = find( trl == trlID(jt) );
%             
%             x = ts(trlIx);
%             
%             [n(jt,:),~] = hist(x,dt);
%             [n2,~] = hist(x,dt2);
%             [xc(jt,:),lag] = xcorr(n2,200);
%             
%             y = jt*ones(1,length(x));
%             x = [x;x];
%             y = [y-.5;y+.5];
%             
%             line(x,y,'Color','k');            
%         end;
%         ylim([0 max(trlID)+1]);
%         xlim([-1e3 7e3]);
%         plot([0 0],[0 max(trlID)+1],'r');
%         plot([2e3 2e3],[0 max(trlID)+1],'r');
%         title(chanLab{BFid(it)});
%         
%         subplot(6,4,[3 4 7 8]);
%         bar(0:200,isi);
%         axis tight;
%         
%         subplot(6,4,[11 12 15 16]);
%         plot(lag((length(lag)-1)/2+2:end),mean(xc(:,(length(lag)-1)/2+2:end),1));
%         axis tight;
%         
%         subplot(6,4,[21 22]);
%         hold on;
%         M = sum(n,1)./size(n,1)./0.25;
%         SE = 1/size(n,1)*sum((n-ones(size(n,1),1)*M).^2,1);
%         SE = SE/sqrt(size(n,1)-1);
%         for kt = 1:length(dt)
%             plot(dt(kt)*ones(1,2),[M(kt) M(kt)+SE(kt)],'k');
%             plot([dt(kt)-50 dt(kt)+50],M(kt)+SE(kt)*ones(1,2),'k');
%             plot([dt(kt)-50 dt(kt)+50],M(kt)+SE(kt)*ones(1,2),'k');
%         end;
%         plot(dt,M,'ks-','MarkerFaceColor','k');
%         xlim([-1e3 7e3]);
%         
%     end;
%     
% end;
% 
% %%
% SUAsel = {};
% MUAsel = {};
% c1 = 0;c2 = 0;
% for it = 1:length(uSel)
%     
%     for jt = 1:length(uSel{it})
%         if ~isempty(uSel{it}) && (uSel{it}(jt)==1)
%             c1 = c1+1;
%             SUAsel{c1} = [it jt];
%         elseif ~isempty(uSel{it}) && (uSel{it}(jt)==2)
%             c2 = c2+1;
%             MUAsel{c2} = [it jt];
%         end;
%     end;
%     
% end;
% 
% x = [MUAsel{:}];
% elecIx = x(1:2:end);
% cluIx = x(2:2:end);
% 
% elecID = unique(elecIx);
% MUAdat = struct;
% for it = 1:length( elecID )
%     
%     [sel] = find(elecIx == elecID(it));
%     
%     ts = [];
%     trl = [];
%     for jt = 1:length( sel )
%         ts = [ts sortedSpikes{sel(jt)}.SpikeTimesSeg];
%         trl = [trl sortedSpikes{sel(jt)}.trl];
%     end;
%     MUAdat.ts{it} = ts;
%     MUAdat.trl{it} = trl;
%     
% end;
% 
% %%
% m = length( MUAdat.trl );
% C = zeros(m,3);
% cIx = 0:255;
% sel = randperm(length(cIx));
% for it = 1:m
%     for jt = 1:3
%         rIx = randperm(length(sel));
%         rIx = rIx(1);
%         C(it,jt) = cIx(sel(rIx));
%         sel(rIx) = [];
%     end;
%     
% end;
% C = C./255;
% 
% dt = -500:4000;
% lag = 500;
% 
% xc = [];
% for kt = 1:52
%     
%     %figure;
%     %hold on;
%     for it = 1:m
%         
%         sel = find( MUAdat.trl{it} == kt );
%         x = MUAdat.ts{it}(sel).*1e3;
%         [n,~] = hist(x,dt);
%         xc(kt,it,:) = xcorr(n,lag);
%         xc(kt,it,lag:lag+2) = NaN;
%         %y = it*ones(1,length(x));
%         %x= [x;x];
%         %y= [y-.5;y+.5];
%         
%         %line(x,y,'Color',C(it,:));
%     end;
%     %xlim([dt(1) dt(end)]);
%     %ylim([0 it+1]);
%     %set(gca,'YTick',[0:it+1]);
%     
%     %figure;
%     %for it = 1:size(xc,1)
%     %    subplot(5,3,it);
%     %    plot(-lag:lag,xc(it,:),'Color',C(it,:));
%     %    axis tight;
%     %    title(it);
%     %end;
%     
% end;
% 
% for it = 1:size(xc,2);figure;bar(-lag:lag,squeeze(mean( xc(:,it,:),1) ), 'k' );xlim([-80 80]);end;
% 
% %%
% 
% ntrl = 52;
% c = 0;
% n = [];
% dt = -500:4e3;
% for it = 1:length(sel)
%     
%     ix = find(uSel{sel(it)});
%     
%     dum = sortedSpikes{sel(it)};
%     
%     for jt = 1:length(ix)
%         
%         c = c+1;
%         
%         ix1 = find(dum.assignedCluster == ix(jt));
%         ts = dum.SpikeTimesSeg.*1e3;
%         
%         for kt = 1:ntrl
%             ix2 = find(dum.trl == kt);
%             ts2 = ts(ix2);
%             ts2(ts2<dt(1)) = [];
%             ts2(ts2>dt(end)) = [];
%             [n(c,kt,:),~] = hist(ts2,dt);
%         end;
%         
%     end;
%     
% end;
% 
% %%
% ix = [ixL ixR];
% 
% [b,a] = butter(4,[3/(Fs/2) 10/(Fs/2)],'bandpass');
% 
% lag = 500;
% dt = -500:1:4000;
% dt2 = -500:4000;
% dt3 = 2000:4000;
% N1 = zeros(1,length( dt ));
% N2 = zeros(1,length( dt ));
% XC1 = [];%zeros(length(-lag:lag),length( dt2 )-lag);
% XC2 = [];%zeros(length(-lag:lag),length( dt2 )-lag);
% Y1 = [];
% Y2 = [];
% clab = [];
% dum2 = [];
% cnt = 0;
% n1 = {};
% n2 = {};
% for kt =1:length( trl )
%     
%     c = 0;
%     ts = {};
%     chck = [];
%     sel = {};
%     for it = 1:length( ix )
%         cID = unique(sortedSpikes{it}.assignedClusterSeg);
%         ix1 = [];
%         ix1 = find( sortedSpikes{it}.trl == trl(kt) );
%         for jt = 1:length( cID )
%             clab = [clab it];
%             if uSel{it}(jt) ==1
%                 
%                 c = c+1;
%                 
%                 sel{c} = it;
%                 
%                 ix2 = [];
%                 ix2 = find( sortedSpikes{it}.assignedClusterSeg == cID(jt) );
%                 
%                 ix3 = intersect(ix1,ix2);
%                 ts{c} = sortedSpikes{it}.SpikeTimesSeg( ix3 ).*1e3;
%                 
%                 if ismember(ix(it),ixL)
%                     chck(c) = 1;
%                 else
%                     chck(c) = 0;
%                 end;
%             end;
%             
%         end;
%         
%     end;
%     
%     %     x = [ts{chck==0}];
%     %     x(x<dt2(1)) = [];
%     %     x(x>dt2(end)) = [];
%     %     nx = hist(x,dt2);
%     %
%     %     x = [ts{chck==0}];
%     %     x(x<dt3(1)) = [];
%     %     x(x>dt3(end)) = [];
%     %     nx2 = hist(x,dt3);
%     %
%     %     Y1(kt,:) = xcorr(nx2,nx2,lag);
%     %
%     %     c = 0;
%     %     xcP1 = [];
%     %     ixT = 1:lag;
%     %     while ixT(end) < length(nx)
%     %         c = c+1;
%     %         [xcP1(:,c)] = xcorr(nx(ixT),nx(ixT),lag);
%     %         %xcP1(500:502,c) = NaN;
%     %         ixT = ixT+1;
%     %     end;
%     %     XC1(kt,:,:) = xcP1;
%     %
%     %     x = [ts{chck==1}];
%     %     x(x<dt2(1)) = [];
%     %     x(x>dt2(end)) = [];
%     %     nx = hist(x,dt2);
%     %
%     %     x = [ts{chck==1}];
%     %     x(x<dt3(1)) = [];
%     %     x(x>dt3(end)) = [];
%     %     nx2 = hist(x,dt3);
%     %
%     %     Y2(kt,:) = xcorr(nx2,nx2,lag);
%     %
%     %     c = 0;
%     %     xcP2 = [];
%     %     ixT = 1:lag;
%     %     while ixT(end) < length(nx)
%     %         c = c+1;
%     %         [xcP2(:,c)] = xcorr(nx(ixT),nx(ixT),lag);
%     %         ixT = ixT+1;
%     %     end;
%     %     XC2(kt,:,:) = xcP2;
%     %
%     %     figure;
%     %     subplot(211);
%     %     imagesc(-250:3750,-lag:lag,xcP1);
%     %     caxis([0 7]);
%     %     subplot(212);
%     %     imagesc(-250:3750,-lag:lag,xcP2);
%     %     caxis([0 7]);
%     
%     sel = [sel{:}];
%     elecID = unique( sel );
%     C1 = [];
%     c = 0;
%     for it = 45:10:255
%         c = c+1;
%         C1(c,:) = [0 0 it]./255;
%     end;
%     rIx1 = randperm(size(C1,1));
%     
%     C2 = [];
%     c=0;
%     for it = 45:10:255
%         c = c+1;
%         C2(c,:) = [it 0 0]./255;
%     end;
%     rIx2 = randperm(size(C2,1));
%     
%     figure;
%     subplot(5,1,1:3);
%     hold on;
%     n1{kt} = [];
%     n2{kt}= [];
%     c1=0;c2 = 0;
%     for it = 1:length( ts )
%         x = [];y = [];
%         x = [ts{it}];
%         y = it*ones(1,length( x ));
%         x = [x;x];
%         y = [y-.5;y+.5];
%         if chck(it) ==1
%             line(x,y,'Color',C2(rIx1(1),:));
%             rIx1(1) = [];
%             dum = ts{it};
%             dum(dum<dt(1)) = [];
%             dum(dum>dt(end)) = [];
%             c1 = c1+1;
%             [n1{kt}(c1,:),~] = hist(dum,dt);
%             
%         else
%             line(x,y,'Color',C2(rIx2(1),:));
%             rIx2(1) = [];
%             dum = ts{it};
%             dum(dum<dt(1)) = [];
%             dum(dum>dt(end)) = [];
%             c2 = c2+1;
%             [n2{kt}(c2,:),~] = hist(dum,dt);
%         end;
%     end;
%     N1 = N1+sum(n1{kt},1);
%     N2 = N2+sum(n2{kt},1);
%     
%     xlim([-500 4000]);
%     subplot(5,1,4);
%     hold on;
%     plot(dt,sum(n1{kt},1),'r-');
%     plot(dt,sum(n2{kt},1),'b-');
%     xlim([-500 4000]);
%     %     subplot(5,1,5);
%     %     hold on;
%     %     for it = 1:length(LFPavg)
%     %         y = filtfilt(b,a,LFPavg{it}(:,kt));
%     %         plot(trlTime.*1e3, y,'b');
%     %     end;
%     %     axis tight;xlim([-500 4000]);
%     
%     %     sel1 = find(chck ==1);
%     %     c = 0;
%     %     ix3 = [];
%     %     for it = 1:length(sel1)
%     %         for jt = it+1:length(sel1)
%     %             c = c+1;
%     %             ix3(c,:) = [it jt];
%     %         end;
%     %     end;
%     %
%     %     dt = -1000:5000;
%     %     lag = 300;
%     %     XC = zeros(length(dt),2*lag+1);
%     %     for it = 1:size(ix3,1)
%     %         ts1 = ts{ix3(it,1)};
%     %         ts2 = ts{ix3(it,2)};
%     %         ts1(logical([(ts1<dt(1))+(ts1>dt(end))]))=[];
%     %         ts2(logical([(ts2<dt(1))+(ts2>dt(end))]))=[];
%     %
%     %         [n1,~] = hist(ts1,dt);
%     %         [n2,~] = hist(ts2,dt);
%     %         xc = zeros(length(dt),2*lag+1);
%     %         parfor jt = 1:length(dt)
%     %             if (jt > lag/2) && (jt < (length(dt)-lag/2))
%     %                 x1 = n1(jt-lag/2:jt+lag/2);
%     %                 x2 = n2(jt-lag/2:jt+lag/2);
%     %                 xc(jt,:) = xcorr(x1,x2,lag);
%     %             end;
%     %         end;
%     %         XC = XC+xc;
%     %     end;
%     %
%     %     sel2 = find(chck ==0);
%     %     c = 0;
%     %     ix2 = [];
%     %     for it = 1:length(sel2)
%     %         for jt = it+1:length(sel2)
%     %             c = c+1;
%     %             ix2(c,:) = [it jt];
%     %         end;
%     %     end;
%     %
%     %     dt = -1000:5000;
%     %     lag = 300;
%     %     XC2 = zeros(length(dt),2*lag+1);
%     %     for it = 1:size(ix2,1)
%     %         ts1 = ts{ix2(it,1)};
%     %         ts2 = ts{ix2(it,2)};
%     %         ts1(logical([(ts1<dt(1))+(ts1>dt(end))]))=[];
%     %         ts2(logical([(ts2<dt(1))+(ts2>dt(end))]))=[];
%     %
%     %         [n1,~] = hist(ts1,dt);
%     %         [n2,~] = hist(ts2,dt);
%     %         xc = zeros(length(dt),2*lag+1);
%     %         parfor jt = 1:length(dt)
%     %             if (jt > lag/2) && (jt < (length(dt)-lag/2))
%     %                 x1 = n1(jt-lag/2:jt+lag/2);
%     %                 x2 = n2(jt-lag/2:jt+lag/2);
%     %                 xc(jt,:) = xcorr(x1,x2,lag);
%     %             end;
%     %         end;
%     %         XC2 = XC2+xc;
%     %     end;
%     %
%     %     figure;
%     %     subplot(121);
%     %     imagesc(dt,-lag:lag,XC');
%     %     xlim([-500 4000]);
%     %     subplot(122);
%     %     imagesc(dt,-lag:lag,XC2');
%     %     xlim([-500 4000]);
%     
% end;
% 
% %%
% AVG = [];
% for zt = 1:length( n1 )
%     x = n1{zt};
%     tsIx = find( sum(x,1) >2 );
%     
%     eIx = {};
%     for it = 1:length( tsIx )
%         
%         eIx{it} = find( x(:,tsIx(it)) ~= 0);
%         
%     end;
%     
%     avg = [];
%     for it = 1:length( eIx )
%         
%         if tsIx(it)-300 >0 & tsIx(it)+300 < size(x,2)
%             dum = x(eIx{it},tsIx(it)-300:tsIx(it)+300);
%             
%             c = 0;
%             eXC = [];
%             for jt = 1:size(dum,1)
%                 x1 = dum(jt,:);
%                 for kt = jt+1:size(dum,1)
%                     
%                     c=c+1;
%                     x2 = dum(kt,:);
%                     eXC(c,:) = xcorr(x1,x2,300);
%                 end;
%             end;
%             avg(it,:) = mean(eXC,1);
%         end;
%     end;
%     
%     if ~isempty( avg )
%         figure;
%         for jt = 1:size( avg,1)
%             subplot(5,4,jt)
%             plot(-300:300,avg(jt,:));
%             axis tight;
%         end;
%     end;
%     
% end;
% 
% 
% 
% %%
% ix1 = find(dt > 0 & dt <= 2e3);
% ix2 = find(dt > 2e3 & dt <= 4e3);
% 
% E = [];
% C = [];
% for it = 1:length( n2 )
%     
%     C(it) = length(find(sum(n1{it}(:,ix1),1)>10));
%     E(it) = length(find(sum(n1{it}(:,ix2),1)>10));
%     
% end;
% 
% figure;
% bar([1 2],[mean(C) mean(E)]);
% 
% %%
% T = size(Y1,2)/1e3;
% W = 1/T;
% TW = T*W;
% k = 2*TW-1;
% params                  = [];
% params.pad              = 2;
% params.Fs               = 1e3;
% params. tapers          = [TW k];
% params.fpass            = [0 30];
% params.trialave         = 1;
% 
% [S,f,R] = mtspectrumpb(Y2',params);

