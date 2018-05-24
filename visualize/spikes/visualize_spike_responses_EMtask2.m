%% include in path
addpath(genpath('/media/rouxf/rds-share/Fred/code/mcode/custom/'));
addpath(genpath('/media/rouxf/rds-share/Fred/code/mcode/toolboxes/wave_clus-testing/'));

%% IDs of patients to include
pID = {'P02' 'P04' 'P05' 'P22AMS' 'P23AMS'};%patient ID

%% loop over patients
for pt = 1%1:length(pID)
    
    fprintf([num2str(pt),'/',num2str( length(pID) )],'\n');
    
    %% set the path for data access
    rpath = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',pID{pt},'/'];%rooth path
    
    [exp,~] = get_EM_rec_info(rpath);% extract session labels
    
    [rpath] = [rpath,exp,filesep]; % path to read data
    
    %%
    pooledDat = dir([rpath,pID{pt},'_pooledSPKdataEMtask.mat']);
    load([rpath,pooledDat.name]);
    
    %%    
    nClu = length([dat.spikes{:}]);
    
    [par]                       = set_parameters_Bham(32000);
    
    colors = [[0.0 0.0 1.0];[1.0 0.0 0.0];[0.0 0.5 0.0];[0.620690 0.0 0.0];[0.413793 0.0 0.758621];[0.965517 0.517241 0.034483];
        [0.448276 0.379310 0.241379];[1.0 0.103448 0.724138];[0.545 0.545 0.545];[0.586207 0.827586 0.310345];
        [0.965517 0.620690 0.862069];[0.620690 0.758621 1.]];
    
    chck=[];for lt = 1:length(dat.sel_b);chck(lt)= isempty(dat.sel_b{lt});end;
    selIdx = find(chck==0);
    
    dat.sel_a = dat.sel_a(selIdx)
    dat.sel_b = dat.sel_b(selIdx);
        
%     figure;
%     c = 0;
%     ax1 = [];
%     for jt = 1:length( dat.spikes )
%         for kt = 1:length( dat.spikes{jt} )
%             c =c+1;
%             subplot(ceil(nClu/4),4,c);
%             ax1(c) = gca;
%         end;
%     end;
%     
%     figure;
%     c = 0;
%     ax2 = [];
%     for jt = 1:length( dat.spikes )
%         for kt = 1:length( dat.spikes{jt} )
%             c =c+1;
%             subplot(ceil(nClu/4),4,c);
%             ax2(c) = gca;
%         end;
%     end;
    
    c = 0;
    mx = [];
    for jt = 1:length( dat.spikes )
        classID = dat.sel_b{jt}+1;
        for kt = 1:length( dat.spikes{jt} )
                        
            maxc = size(colors,1);
            class = dat.classes{jt}{kt};
            nclasses = max(class);
                        
            [~, SCORE] = pca( dat.spikes{jt}{kt} );
            aux_spk = find( class == classID(kt) );
            
            
            c = c+1;
            %axes(ax1(c));            
            %subplot(ceil(nClu/4),4,c);
            figure;
            hold on;
            plot(SCORE(:,1),SCORE(:,2),'.','Color',[.5 .5 .5],'MarkerSize',.5);                        
            plot(SCORE(aux_spk,1),SCORE(aux_spk,2),'.','Color',colors(mod(c-1,maxc)+1,:),'MarkerSize',.5);
            axis tight;axis off;
            set(gcf,'Color','w');
            
            %axes(ax2(c));
            %subplot(ceil(nClu/4),4,c);   
            figure;
            hold on;
            plot(1/32000:1/32000:0.002,dat.spikes{jt}{kt}(aux_spk,:),'Color',colors(mod(c-1,maxc)+1,:));
            plot(1/32000:1/32000:0.002,mean(dat.spikes{jt}{kt}(aux_spk,:),1),'Color',[0 0 0],'LineWidth',3);
            axis tight;axis off;
            
            mx(c,1) = [min(min(dat.spikes{jt}{kt}(aux_spk,:)))];
            mx(c,2) = [max(max(dat.spikes{jt}{kt}(aux_spk,:)))];
            
            plot([0.004 0.004],[mx(c,1) mx(c,1)+100],'k','LineWidth',3);
            plot([0.003 0.004],[mx(c,1) mx(c,1)],'k','LineWidth',3);
            set(gcf,'Color','w');
            
        end;
    end;
    %set(ax2,'YLim',[min(min(mx)) max(max(mx))]);
    
    %%
    % [inspk]                     = wave_features(spikes,par);
    % inputs = size(inspk,2);
    %
    % figure(10)
    % for i=1:inputs
    %     for j=i+1:inputs
    %         %subplottight(inputs,inputs,(i-1)*inputs+j)
    %         %subplot(inputs-1,inputs-1,(i-1)*(inputs-1)+j-1)
    %         %figure;
    %         [c,r] = ind2sub([inputs inputs], (i-1)*inputs+j);
    %         ax = subplot('Position', [(c-1)/inputs, 1-(r)/inputs, 1/inputs, 1/inputs]);
    %         hold on
    %         for k=1:nclasses
    %             class_aux = find(classes==k);
    %             max_spikes = min(par.max_spk,length(class_aux));
    %             plot(inspk(class_aux(1:max_spikes),i),inspk(class_aux(1:max_spikes),j),'.','color', colors(mod(k-1,maxc)+1,:),'markersize',.5)
    %             axis off
    %         end
    %     end
    % end
    
    %% make raster plot
    for xt = 1:length( dat.cluRaster )
        for jt = 1:length( dat.cluRaster{xt} )
            x = dat.cluRaster{xt}{jt};
            figure;
            subplot(5,1,1:3);
            hold on;
            plot([0 0],[0 size(x,1)+1],'r');
            plot([2000 2000],[0 size(x,1)+1],'r');
            t = {};
            for it = 1:size(x,1)
                
                ix = find(x(it,:)~= 0);
                if  ~isempty(ix)
                    t{it} = dat.dt(ix);
                    d = [t{it};t{it}];
                    d2 = it*ones(1,length(ix));
                    d2 = [d2-.5;d2+.5];
                    
                    line(d,d2,'Color','k');
                end;
            end;
            xlim([dat.dt(1) dat.dt(end)]);
            set(gcf,'Color','w');
            set(gca,'YLim',[0 length(dat.trl_all{xt})+1]);
            title(dat.CSC{xt}(dat.sel_a{xt}(jt)));
            
            dt2 = dat.dt(1):500:dat.dt(end);
            n = zeros( length(t) , length(dt2) );
            for it = 1:length( t )
                if ~isempty(t{it})
                    n(it,:) = hist( t{it} , dt2 );
                else
                    n(it,:) = zeros( 1 , length( dt2 ) );
                end;
            end;
            n = n./(500/dat.Fs);
            n(:,[1 end]) = [];
            dt2([1 end]) = [];
            
            M = sum(n,1)./size(n,1);
            SD = std(n,0,1)./sqrt(size(n,1)-1);
            
            subplot(5,1,4:5);
            hold on;
            plot([0 0],[min(M-SD) max(M+SD)],'r');
            plot([2000 2000],[min(M-SD) max(M+SD)],'r');
            errorbar(dt2+250,M,SD,'ks-','MarkerFaceColor','k');
            axis tight;xlim([dat.dt(1) dat.dt(end)]);
            set(gcf,'Color','w');
       end;
        
    end;
    
%     %% plot mean firing rate              
%     figure;
%     subplot(221);
%     a = gca;
%     hold on;
%     cond = [1 2];
%     for it = 1:length(dat.FRc)
%         for jt = 1:length( cond )
%             plot(jt,dat.FRc{it}(cond(jt)),'o','Color','b','MarkerFaceColor','k');
%         end;
%     end;
%     for it = 1:length( dat.FRc )
%         plot([1 2],dat.FRc{it}(cond),'k-');
%     end;
%     
%     subplot(223);
%     a = [a gca];
%     hold on;
%     cond = [1 2];
%     for it = 1:length(dat.FRe)
%         for jt = 1:length( cond )
%             plot(jt,dat.FRe{it}(cond(jt)),'o','Color','r','MarkerFaceColor','k');
%         end;
%     end;
%     for it = 1:length( dat.FRe )
%         plot([1 2],dat.FRe{it}(cond),'k-');
%     end;
%     
%     subplot(222);
%     a = [a gca];
%     hold on;
%     n = zeros(length(dat.selIdxC),2);
%     for it = 1:size(dat.nItemsC,1)
%         for jt = 1:size(dat.nItemsC,2)
%             plot(jt,dat.nItemsC(it,jt),'o','Color','b','MarkerFaceColor','k');
%         end;
%     end;
%     for it = 1:size(dat.nItemsC,1)
%         plot([1 2],dat.nItemsC(it,:),'k-','MarkerFaceColor','k');
%     end;
%     
%     subplot(224);
%     a = [a gca];
%     hold on;
%     n = zeros(length(dat.selIdxE),2);
%     for it = 1:size(dat.nItemsE,1)
%         for jt = 1:size(dat.nItemsE,2)
%             plot(jt,dat.nItemsE(it,jt),'o','Color','b','MarkerFaceColor','k');
%         end;
%     end;
%     for it = 1:size(dat.nItemsE,1)
%         plot([1 2],dat.nItemsE(it,:),'k-','MarkerFaceColor','k');
%     end;
%     
%     for it = 1:length(a)
%         box(a(it),'on');
%     end;
%     set(a,'Xlim',[0 3]);
%     set(a(1),'XTickLabel',{'' 'Baseline' 'Cue' ''});
%     set(a(2),'XTickLabel',{'' 'Baseline' 'Encoding' ''});
%     set(a(3:4),'XTickLabel',{'' 'Low' 'High' ''});
%     set(gcf,'Color','w');
%     ylabel(a(1),'Firing rate (Hz)');
%     ylabel(a(2),'Firing rate (Hz)');
%     ylabel(a(3),'Mean item recall');
%     ylabel(a(4),'Mean item recall');
    
end;



%%
[spkTms1,cond1] = clusterSpikeTimes2mat(sel_ix1,sel_ix2,spike_dat,trl_all{xt},ix(1:3),[-2000 5000],250);
[spkTms2,cond2] = clusterSpikeTimes2mat(sel_ix3,sel_ix4,spike_dat,trl_all{xt},ix(1:3),[-2000 5000],250);
[spkTms3,cond3] = clusterSpikeTimes2mat(sel_ix5,sel_ix6,spike_dat,trl_all{xt},ix(1:3),[-2000 5000],250);
[spkTms4,cond4] = clusterSpikeTimes2mat(sel_ix3,sel_ix4,spike_dat,trl_all{xt},ix(4:6),[-2000 5000],250);

%%
visualize_SUrespEM(spike_dat,sel_ix1,sel_ix2,spkTms1,cond1,nLog,trl_all{xt},-2000:250:5000);
visualize_SUrespEM(spike_dat,sel_ix3,sel_ix4,spkTms2,cond2,nLog,trl_all{xt},-2000:250:5000);

%%
[h] = visualizeFRcondAVG(sel_ix3,spkTms2,cond1,-2000:250:5000,1000,[2000 5000]);
[h] = visualizeFRcondAVG(sel_ix3(3),spkTms4(3),cond4,-2000:1:5000,1000,[2000 5000]);

%% waveform plots
color = {'c' 'm' 'g' 'y' 'r' 'b'};
c = 0;
for it = 1:length(sel_ix3)
    c =c+1;
    x = squeeze( spike_dat{sel_ix3(it)}.waveform{sel_ix4(it)+1} );
    if c > length(color)
        c = 1;
    end;
    
    figure;
    hold on;
    plot(spike_dat{sel_ix3(it)}.waveformtime, x , 'Color',color{c});
    plot(spike_dat{sel_ix3(it)}.waveformtime, mean(x,2) , 'k','LineWidth',3);
    axis tight;
    %plot([2e-3 2e-3],[min(min(x)) min(min(x))+25],'k','LineWidth',3);
    %plot([1.5e-3 2e-3],[min(min(x)) min(min(x))],'k','LineWidth',3);
    box off;axis off;
    set(gcf,'Color','w');
    
    %     ts = spike_dat{sel_ix3(it)}.timestamp{sel_ix4(it)+1};
    %     ts = ts./1e3;
    %     isi = diff(ts);
    %     dt = 0:5:500;
    %     n = histc(isi,dt);
    
    %     figure;
    %     bar(dt,n);
    %     axis tight;
    %     xlabel('ISI [ms]');
    %     ylabel('Count');
    %     box off;
    %     set(gcf,'Color','w');
    
end;

%%
[sel_ix] = find(CSC(:,end)==1);
i_base = find( dt >= -1000 & dt < 0 );
i_cue1 = find( dt >= 0 & dt < 1000);
i_cue2 = find( dt >= 1000 & dt < 2000);
i_enc1 = find( dt >= 2000 & dt < 3000);
i_enc2 = find( dt >= 3000 & dt < 4000);

for jt = 1:length(sel_ix)
    
    x = n{sel_ix(jt)};
    
    sc1 = sum(x(i_base,:),1);
    sc2 = sum(x(i_cue1,:),1);
    sc3 = sum(x(i_cue2,:),1);
    sc4 = sum(x(i_enc1,:),1);
    sc5 = sum(x(i_enc2,:),1);
    
    figure;
    hold on;
    for it = 1:5;
        h(it) = bar(it,eval(['mean(sc',num2str(it),')']));
    end;
    
    set(h(1),'FaceColor','k')
    set(h(2:3),'FaceColor','r')
    set(h(4:5),'FaceColor','b')
    set(h(1),'FaceColor',[.75 .75 .75])
    set(h,'LineWidth',3);
    
    for it = 1:5
        plot([it it],[mean(eval(['sc',num2str(it)]))-std(eval(['sc',num2str(it)]))/sqrt(length(eval(['sc',num2str(it)]))-1) mean(eval(['sc',num2str(it)]))+std(eval(['sc',num2str(it)]))/sqrt(length(eval(['sc',num2str(it)]))-1)],'k','LineWidth',3)
        plot([it-.1 it+.1],[mean(eval(['sc',num2str(it)]))-std(eval(['sc',num2str(it)]))/sqrt(length(eval(['sc',num2str(it)]))-1) mean(eval(['sc',num2str(it)]))-std(eval(['sc',num2str(it)]))/sqrt(length(eval(['sc',num2str(it)]))-1)],'k','LineWidth',3)
        plot([it-.1 it+.1],[mean(eval(['sc',num2str(it)]))+std(eval(['sc',num2str(it)]))/sqrt(length(eval(['sc',num2str(it)]))-1) mean(eval(['sc',num2str(it)]))+std(eval(['sc',num2str(it)]))/sqrt(length(eval(['sc',num2str(it)]))-1)],'k','LineWidth',3)
        
    end;
    set(gca,'XTick',[1 2.5 4.5]);
    set(gca,'XTickLabel',{'Baseline' 'Cue' 'Encoding'});
end;

%%
mwID = unique(CSC(:,2));

nc = [];
for it = 1:length(mwID)
    
    ix = find(CSC(:,2) == mwID(it));
    nc(it) = length(ix);
    
end;


[nc,xAx] = hist(nc);

figure;bar(xAx,nc);

%%
sel_ix = find(CSC(:,end)==0);

AVG1 = zeros(size(Cgrm{1}));
for it = 1:length(sel_ix)
    
    Cgrm{sel_ix(it)}(isnan(Cgrm{sel_ix(it)})) = 0;
    
    AVG1 = AVG1 + Cgrm{sel_ix(it)};
    
end;
AVG1 = AVG1./length(sel_ix);

sel_ix = find(CSC(:,end)==1);
AVG2 = zeros(size(Cgrm{1}));
for it = 1:length(sel_ix)
    
    AVG2 = AVG2 + Cgrm{sel_ix(it)};
    
end;
AVG2 = AVG2./length(sel_ix);

figure;
subplot(211);
hold on;
imagesc(tAx-2,fAx,AVG1');
axis xy;
plot([0 0],[min(fAx) max(fAx)],'w');
plot([2 2],[min(fAx) max(fAx)],'w');
axis tight;

subplot(212);
hold on;
imagesc(tAx-2,fAx,AVG2');
axis xy;
plot([0 0],[min(fAx) max(fAx)],'w');
plot([2 2],[min(fAx) max(fAx)],'w');
axis tight;