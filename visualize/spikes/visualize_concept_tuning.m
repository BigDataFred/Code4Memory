%% analysis parameters 

basepath = '/media/rouxf/rds-share/iEEG_DATA/MICRO/';
%basepath = '/media/rouxf/rds-share/iEEG_DATA/MICRO/tmp/';
%basepath = '/home/rouxf/out/';% path from which to read data

pID = 'P09';% patient ID
refP = [0.003 3];% pct for SU ISI check
mode = 'sorted';% use sorted vs unsorted data
algo = 'waveclus';% select between 'waveclus' and 'osort'

elec_sel = {};%{'antHippL7'};%{'midHippR6','midHippR8'};% % examine elec of interest
cluster_sel =  [];%[2 2];%[2];%[ 2 ];% examine cluster of interest

%% extract the session labels from directory

d = dir([basepath,pID,filesep,'Tunings',filesep]);
%d = dir([basepath,pID,filesep]); % 

d(1:2) = [];

sesh_ts = cell(1,length(d));
for it = 1:length(d)
    sesh_ts(it) = {d(it).name};
end;
sesh_ts = sesh_ts'

%% loop over the different session
for sesh_it = 4%1:length(sesh_ts)
    
    sesh_it
    
%     if sesh_it <6
        bt = [ 0.55 0.055 ];%500ms
        pt = [ 0.2 0.7 ];%500ms
%     else
%         bt = [ 0.9 0.1 ];%600ms
%         pt = [ 0.3 1.2 ];%600ms
%     end;
    
    p2sf = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',pID,'/Tunings/',sesh_ts{sesh_it},'/spike_dat/'];
    p2lf = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',pID,'/Tunings/',sesh_ts{sesh_it},'/log_dat/'];
    
    %p2sf = ['/home/rouxf/out/',pID,'/',sesh_ts{sesh_it},'/spike_dat/'];
    %p2lf = ['/home/rouxf/out/',pID,'/',sesh_ts{sesh_it},'/log_dat/'];
    
    CSC_files = dir([p2sf,pID,'_*_spike_data_CSC_*_stimlocked.mat']);
    LOG_file  = dir([p2lf,pID,'_*_log_ctune_*_LogDat.mat']);
    load([p2lf,LOG_file(end).name]);
    
    %% in case an electrode of interest has been specified
    if ~isempty(elec_sel)
        selCSCix = [];
        for zt = 1:length(elec_sel)
            chck = regexp({CSC_files(:).name}',elec_sel{zt});
            for yt = 1:length( chck )
                if ~isempty(chck{yt} )
                    selCSCix(zt) = yt;
                end;
            end;
        end;
        CSC_files = CSC_files( selCSCix );
    end;
    
    %% make an exception for patients P04 and P05
    if strcmp(pID,'P04') || strcmp(pID,'P05')
        
        ID1 = zeros(length(CSC_files)/2,2);
        ID2 = zeros(length(CSC_files)/2,2);
        
        c1 = 0; c2 = 0;
        for jt = 1:length(CSC_files)
            
            ix = regexp(CSC_files(jt).name,'CSC_R\w*\d')+5:regexp(CSC_files(jt).name,'_2017')-1;
            
            for it = 1:24
                
                if strcmp(['CSC_R',num2str(CSC_files(jt).name(ix))],['CSC_R',num2str(it)])
                    c1 = c1+1;
                    ID1(c1,1) = jt;
                    ID1(c1,2) = str2double( CSC_files(jt).name(ix) )+24;
                end;
            end;
            
            ix = regexp(CSC_files(jt).name,'CSC_L')+5:regexp(CSC_files(jt).name,'_2017')-1;
            
            for it = 1:24
                if strcmp(['CSC_L',num2str(CSC_files(jt).name(ix))],['CSC_L',num2str(it)])
                    c2 = c2+1;
                    ID2(c2,1) = jt;
                    ID2(c2,2) = str2double( CSC_files(jt).name(ix) );
                end;
            end;
            
        end;
        
        del_ix = find(ID1(:,1) ==0);
        ID1(del_ix,:) = [];
        
        del_ix = find(ID2(:,1) ==0);
        ID2(del_ix,:) = [];
        
        sel_ix = [ID2(:,1);ID1(:,1)]
        
        if ~isequal({CSC_files(sel_ix).name}',{CSC_files(:).name}')
            error('index labels must match');
        end;
        
        ix = [ ID2(:,2)' ID1(:,2)' ];
        
        [~,s_idx] = sort( ix );
        
        ix = [ ID2(:,1)' ID1(:,1)' ];
        
        ix = ix( s_idx )';
        
        CSC_files = CSC_files(ix);
    end;
            
    %% define colors for plotting        
    C = [];
    C(1,:) = [255 255 0];%Y
    C(2,:) = [102 255 255];%C
    C(3,:) = [178 102 255];%M
    C(4,:) = [255 153 51];%O
    C(5,:) = [128 255 0];%G
    C(6,:) = [255 102 102];%R
    C = C./255;        
    
    %% make PCA plot for each BF electrode
    figure;
    CS = {};
    for gt = 1:length(CSC_files)% loop over the BFs
        
        load([p2sf,CSC_files(gt).name]);% load the data of each BF
        
        %%
        
        switch algo
            case 'waveclus'
                switch mode
                    case 'sorted'
                        spike_data = save_data{1}{1}{1};
                        
                    case 'unsorted'
                        spike_data = save_data{1}{1}{2};
                        
                end;
            case 'osort'
                switch mode
                    case 'sorted'
                        spike_data = save_data{1}{1}{2};
                        
                    case 'unsorted'
                        spike_data = save_data{1}{1}{4};
                        
                end;
        end;
        
        %%       
        SCORE = cell(1,length(spike_data.waveform));
        for it = 1:length(SCORE)
            
            [~, SCORE{it}] = pca(squeeze(spike_data.waveform{it})');
            if ~isempty(SCORE{it})
                subplot(length(CSC_files)/8,8,gt);
                hold on;
                if size(SCORE{it},2)>1
                    plot(SCORE{it}(:,1),SCORE{it}(:,2),'.','Color',[.75 .75 .75]);
                end;
            end;
            
        end;
        
        %%
        if ~isempty([spike_data.unit{:}])
            cID   = 1:length(spike_data.unit);
            cID(cID==1) = [];
            
            for jt = 1:length( cID )
                                
                X = spike_data.waveformtime.*1e3;
                trl = spike_data.trial{ cID(jt) };
                
                ts = (double(spike_data.timestamp{ cID(jt) }).*1e-6).*1e3;
                
                dt = [0:1:200 201];
                
                lm = [];
                Y = cell(1,length(cID));
                
                f =0; mem = [];
                thr = {}; thr2 = {};
                c= 0;
                for it = 1:length(cID)
                    c = c+1;
                    Y{it} = squeeze( spike_data.waveform{ cID(jt) } );
                    
                    if ~isempty(Y{it})
                        lm = [lm min(min(Y{it})) max(max(Y{it}))];
                    end;                    
                                                           
                    subplot(length(CSC_files)/8,8,gt)
                    hold on;      
                    if size(SCORE{jt},2)>1
                        plot(SCORE{jt}(:,1),SCORE{jt}(:,2),'.','Color',C(jt,:));                    
                    end;
                                        
                end;
                
                axis tight;
                xlabel('PC 1');
                ylabel('PC 2');
            end;
        end;
    end;
    
    %% extract the picture event data from LogFile
    picev = LogDat.dat(:,2);% file names of each image
    picev = unique(picev);% keep each filename once
    
    % get the trial indexes of each piture event
    picix = cell(1,length(picev));
    nrep = zeros(1,length(picev));
    chck = zeros(1,length(picev));
    for it = 1:length(picev)% loop over each image once
        % store the indexes of the trials at which image was presented
        picix{it} = find(strcmp(LogDat.dat(:,2),picev(it)));
        nrep(it) = length(picix{it});% number of times image was presented        
        % check that all trial indexes correspond with the correct image                
        chck(it) = sum(strcmp( LogDat.dat(picix{it},2),picev(it) ));
    end;
    nrep = unique(nrep);
    
%     %sanity checks
%     if sum(chck==nrep) ~= length(picev)
%         error('event repetitions must be equal across events');
%     end;
    
%     if length(nrep) >1
%         error('event repetition must be unique');
%     end;
        
%     if ( nrep*length(picev) ) ~= size( LogDat.dat,1)
%         error( 'event repetitions must correspond with trial number' );
%     end;
    
    %% intialize some variables
    pct = cell( 1 , length(CSC_files) );
    nspk = cell( 1 , length(CSC_files) );
    tinf = [];
    chck = [];
    
    %% examine tuning responses    
    for ut = 1:length( CSC_files )% loop over the BF electrodes
        ut
        load([p2sf,CSC_files(ut).name]);% load the data
        
        switch algo
            case 'waveclus'
                switch mode
                    case 'sorted'
                        spike_data = save_data{1}{1}{1};
                        
                    case 'unsorted'
                        spike_data = save_data{1}{1}{2};
                        
                end;
            case 'osort'
                switch mode
                    case 'sorted'
                        spike_data = save_data{1}{1}{2};
                        
                    case 'unsorted'
                        spike_data = save_data{1}{1}{4};
                        
                end;
        end;                      
        
        %% in case a cluster of interest has been pre-specified
        if ~isempty(cluster_sel) && ~isempty(elec_sel)
            
            [selIX] = cluster_sel(ut);
            
            spike_data.label = {spike_data.label{selIX}};
            spike_data.timestamp = {spike_data.timestamp{selIX}};
            spike_data.waveform = {spike_data.waveform{selIX}};
            spike_data.unit = {spike_data.unit{selIX}};
            spike_data.time = {spike_data.time{selIX}};
            spike_data.trial = {spike_data.trial{selIX}};
        end;
        
        %% cluster identity label
        cID = [spike_data.unit{:}];%these are the cluster labels for each CSC channel          
         
        %%
        if ~isempty(cID)
            
            cID = unique(cID);% add 1 because of cluster zero
            
            %initialize, tinf stands for tuning information
            n_isi = []; pct = []; FR =[]; sdf =[]; R ={}; SC = []; tinf ={}; x1 = {}; x2 = {};
            
            for it = 1:length(cID)% loop over the clusters
                
                c = 1;%counter
                
                % tuning information, 1st entry = BF index
                tinf{it}(c) = it; % save the cluster index                                
                
                [st] = spike_data.time{ it };% get spike times for each cluster
                
                [trl] = spike_data.trial{ it };% get the trial labels of each cluster
                
                dum = [];for ot =1:length(picix);dum = [dum picix{ot}'];end; dum = sort(dum(:));
                trlID = dum;% get each trial label only once
                
                % 1) define the baseline spiking activity
                
                dt = [-bt(1):1e-3:-bt(2) -bt(2)+1e-3];
                isi = [];
                x1{it} = []; x2{it} = [];
                for gt = 1:length(picix)% loop over the trial labels
                    n = zeros(1,length(picix{gt}));% count the number of spikes in the baseline of each trial
                    n2 = zeros(length(picix{gt}),length(dt));
                    for ft = 1:length(picix{gt})
                        ix = find(trl == picix{gt}(ft));% indexes of the st corresponding to each trl
                        
                        dum = st( ix );
                        if length( dum ) >1
                            isi = [isi diff(dum)];
                        end;
                        
                        base = dum(dum >=-bt(1) & dum <=-bt(2));%get the sts in the baseline
                        [base,~] = sort(base);
                        n(ft) = length( base );% number of sts
                        [n2(ft,:),~] = hist(base,dt);% make hist of sts
                    end;
                    n2(:,end)= [];
                    x1{it}(gt) = median(n);
                    x2{it}(gt) = median(sum(n2,2));
                end;
                
                isi = diff(st.*1e3);
                isi(sign(isi)==-1) = [];
                
                pct(it) = (length(find(isi <=3))/length(isi))*100;                
                
                dt_isi = 0:1:1501;
                n_isi(it,:) = hist(isi,dt_isi);
                
                if pct(it) < refP(2) || (~isempty(elec_sel))
                
                % define the treshold for baseline spiking activity above which
                % poststim must be elevated
                [ thr{it}  ] = median(x1{it}) + (5*std(x1{it}));
                [ thr2{it} ] = median(x2{it}) + (5*std(x2{it}));% mean number of spikes in the baseline + 5*SDs
                [thr{it} thr2{it}]
                
                if ( thr{it} ~= thr2{it} )%just a sanity check
                    error('threshold values must match');
                end;
                
                % compute spiking activity during poststim
                for jt = 1:length(picix)% loop over stimuli
                    
                    chck = zeros(1,4);% reset the detection flags after each pass
                    
                    rep = picix{jt};% indices at which trial which stimulus was shown
                    
                    RT = LogDat.RT(picix{jt});%reaction times
                    
                    %another sanity check
                    if ( length( unique(LogDat.dat(rep,2)) ) > 1 );
                        error( 'picture events must be unique' );
                    end;
                    
                    B = []; B2 = []; PS = []; PS2 = []; fr = [];
                    for kt = 1:length(rep) %loop over repetitions
                        
                        sel = find( trl == rep(kt));
                        x= st( sel ); % spike times per stim
                        
                        %R{it,jt,kt} = x( x >= -bt(1) & x <= pt(2));
                        R{it,jt,kt} = x( x >= -1 & x <= 1);
                        
                        % compute the mean number of spikes in the baseline across all trials
                        base = x(x>=-bt(1) & x <=-bt(2));% get the spike times during the baseline
                        dt = [-bt(1):1e-3:-bt(2) -bt(2)+1e-3];
                        [n,~] = hist( base,dt );
                        n(end) = [];
                        B(kt) = length( base );% absolute number of spikes
                        B2(kt) = sum(n);% sum over time
                        
                        post = x(x>=pt(1) & x <=pt(2));% get the spike times during the baseline
                        dt = [pt(1):1e-3:pt(2) pt(2)+1e-3];
                        [n,~] = hist( post,dt );
                        n(end) = [];
                        PS(kt) = length( post );% absolute number of spikes
                        PS2(kt) = sum(n,2);%sum over time
                        
                        fr(kt,:) = n;
                    end;
                    
                    T = pt(2)-pt(1);
                    FR(it,jt) = median(sum(fr,2)/T);
                    
                    dum = PS;
                    dum(dum==Inf) = 0;
                    %dum = (dum - mean(x2))./mean(x2);
                    SC(it,jt,1) = nanmedian(dum);
                    SC(it,jt,2) = nanstd(dum)/sqrt(length(dum)-1);
                    
                    sn = 0.1*1000;
                    sn = sn/2;
                    m = 0;
                    sd = 1;
                    gw = pdf('norm',-sn:sn,m,sd);
                    dt = [-1:5e-2:1];
                    x = zeros(size(R,3),length(dt));
                    n = zeros(size(R,3),length(dt));
                    for kt = 1:size(R,3)
                        [n(kt,:),~] = hist(R{it,jt,kt},dt);
                        x(kt,:) = conv(n(kt,:),gw,'same');
                    end;
                    
                    sdf(it,jt,:) = mean(x,1);
                    
                    if ~isequal(PS,PS2) || ~isequal(B,B2)
                        error('spike count must match');
                    end;
                    
                    %thr = mean(B)+(5*std(B));
                    
                    %% test significance of tuning response
                    chck = [];
                    
                    % 1) check that the median number of spikes in the poststim period
                    % is higher than the threshold
                    if ( (SC(it,jt,1)-SC(it,jt,2)) > thr2{it} )%( SC(it,jt,1) > thr2{it} )%
                       chck(1) = 1;
                    end;
                    
                    % 2) check that the median number of spikes in the poststim
                    % interval was at least 2
                    if ( median(PS) >= 2 ) % ( FR(it,jt) >= 2 ) %
                        chck(2) = 1;
                    end;
                    
                    % 3) check that the number of spikes in the PS period is
                    % significantly higher than the number of spikes in the B period
                    [~,p] = ttest( B' , PS');%[p,~] = ranksum( B' , PS');%
                    if ( p < 0.05 )
                        chck(3) = 1;
                    end;
                    
                    if ~isempty(elec_sel) %&& ~isempty(cluster_sel)
                        chck(1:3)=1;
                    end;
                    
                    %% extract the BF-label
                    ix = regexp(CSC_files(ut).name,'CSC_\w{1,2}\d{1,2}');
                    chan = CSC_files(ut).name(ix:ix+regexp(CSC_files(ut).name(ix:end),'\d{1,2}_2016'));
                    chan(regexp(chan,'_')) = [];                    
                    if isempty(ix)
                        ix = regexp(CSC_files(ut).name,'CSC')+4;
                        chan = CSC_files(ut).name(ix:ix+regexp(CSC_files(ut).name(ix:end),'\d{1,2}_201'));
                        chan(regexp(chan,'_')) = [];
                    end;
                
                    %% store information for visualization of significant response
                    if ( sum(chck)== 3 )
                        c = c+1;
                        fprintf(['Tuning on channel: ',chan,'\n']);
                        fprintf(['Cluster: ', num2str( cID(it) ),'\n']);
                        fprintf(['Image: ', cell2mat(unique([LogDat.dat(picix{jt},2)])) ,'\n']);
                        fprintf(['Index: ', num2str(find(strcmp(picev,unique([LogDat.dat(picix{jt},2)])))) ,'\n']);
                        tinf{it}(c) = jt;% save the image event index
                    end;
                    
                end;
                end;
            end;
            
            %% check for clusters with significant respones, and mark them 
            sel = [];
            c = 0;
            for jt = 1:length(tinf) % for 
                if length(tinf{jt})>1
                    c = c+1;
                    sel(c) = jt; % save the indexes of clusters with tunings
                end;
            end;
            
            %% do the visualization of significant tuning responses
            if ~isempty(sel)                                
                
                %% PCA plot
                figure;
                n = 2+length(sel);
                subplot(2,n,[1 2 n+1 n+2]);
                ax1 = gca;
                hold(ax1,'on');
                SCORE = cell(1,length(spike_data.waveform));
                for kt = 1:length( SCORE )
                [~, SCORE{kt}] = pca(squeeze(spike_data.waveform{kt})');
                    if size(SCORE{kt},2)>1
                        plot(ax1,SCORE{kt}(:,1),SCORE{kt}(:,2),'.','Color',[.75 .75 .75]);
                    end;
                end;
                xlabel(ax1,'PC1');
                ylabel(ax1,'PC2');
                axis(ax1,'tight');
                set(ax1,'XTick',[min(get(gca,'XTick')) max(get(gca,'XTick'))]);
                set(ax1,'YTick',[min(get(gca,'YTick')) max(get(gca,'YTick'))]);
                title(ax1,chan);
                
                ix = 3:n;
                ix2 = n+3:2*n;
                ax2 = zeros(1,length(ix));
                ax3 = zeros(1,length(ix2));
                for ot = 1:length(ix)
                    subplot(2,n,ix(ot));
                    ax2(ot) = gca;
                    hold(ax2(ot),'on');
                    subplot(2,n,ix2(ot));
                    ax3(ot) = gca;
                    hold(ax3(ot),'on');
                end;
                
                %% initialize some figures
                n = length([tinf{sel(:)}])-length(sel);                
                
                figure;                
                ax6 = zeros(1,length(sel));
                for ot = 1:length(sel)
                    subplot(1,length(sel),ot);
                    ax6(ot) = gca;
                    hold(ax6(ot),'on');
                end;
                
                figure;
                ax7 = zeros(1,n);
                for ot = 1:n
                    subplot(1,n,ot);
                    ax7(ot) = gca;
                    hold(ax7(ot),'on');
                end;
                
                %% make raster plot                
                cnt1 = 0;
                cnt2 = 0;
                for jt = 1:length(sel)
                    
                    for zt = 2:length(tinf{sel(jt)})
                        
                        cnt1 = cnt1+1;       
                        
                        figure
                        subplot(211);
                        hold on;
                        k = 0;
                        k = k+1;
                        for nt = 1:size(R,3)
                            x = R{tinf{sel(jt)}(1),tinf{sel(jt)}(zt),nt};
                            y = nt*ones(1,length(x));
                            x = [x;x];
                            y = [y-.5;y+.5];
                            
                            plot(x,y,'Color','b','LineWidth',3);
                            %plot(ax4(cnt1),x,y,'Color','b','LineWidth',3);
                        end;
                        
                        imgLabel = unique({LogDat.dat{picix{tinf{sel(jt)}(zt)},2}});
                        imgLabel2 = picev{tinf{sel(jt)}(zt)};
                        if ~strcmp(imgLabel,imgLabel2)
                            error('event label assignment must match');
                        end;
                        [~,pev,~] = fileparts( imgLabel{:} );
                        pev(regexp(pev,'_')) = [];                   
                        imgIx = find(strcmp(picev,unique([LogDat.dat(picix{tinf{sel(jt)}(zt)},2)])));
                        imgIx2 = find(strcmp(picev,imgLabel));
                        if imgIx ~=imgIx2
                            error('event indexes must match');
                        end;
                        title([pev,' (',num2str(imgIx2),')']);
                        axis('off');
                        xlim([-.2 1]);%xlim([-bt(1) pt(2)]);
                        ylim([0 size(R,3)+1]);
                        %title(ax4(cnt1),pev);
                        %axis(ax4(cnt1),'off');
                        %xlim(ax4(cnt1),[-bt(1) pt(2)]);
                        %ylim(ax4(cnt1),[0 size(R,3)+1]);
                        
                        subplot(212);
                        hold on;
                        %tx = -bt(1):5e-2:pt(2);
                        tx = -1:5e-2:1;
                        plot([0 0],[0 max(max(sdf(tinf{sel(jt)}(1),tinf{sel(jt)}(zt),:)))],'k--','LineWidth',3);
                        %plot(ax5(cnt1),[0 0],[0 max(max(sdf(tinf{sel(jt)}(1),tinf{sel(jt)}(zt),:)))],'k--','LineWidth',3);
                        %bar(tx,sdf(tinf(sel(jt),3),:));
                        ix = find(tx >=pt(1) & tx <=pt(2));
                        h = area(tx(ix),squeeze(sdf(tinf{sel(jt)}(1),tinf{sel(jt)}(zt),ix)));
                        %h = area(ax5(cnt1),tx(ix),squeeze(sdf(tinf{sel(jt)}(1),tinf{sel(jt)}(zt),ix)));
                        set(h,'FaceColor',[.75 .75 .75],'EdgeColor',[.75 .75 .75]);
                        txt = round(FR(tinf{sel(jt)}(1),tinf{sel(jt)}(zt))*100)/100;
                        ix = find(tx >=-bt(1) & tx <=-bt(2));
                        h = area(tx(ix),squeeze(sdf(tinf{sel(jt)}(1),tinf{sel(jt)}(zt),ix)));
                        %h = area(ax5(cnt1),tx(ix),squeeze(sdf(tinf{sel(jt)}(1),tinf{sel(jt)}(zt),ix)));
                        set(h,'FaceColor',[.75 .75 .75],'EdgeColor',[.75 .75 .75]);
                        plot(tx,squeeze(sdf(tinf{sel(jt)}(1),tinf{sel(jt)}(zt),:)),'r','LineWidth',3);
                        %plot(ax5(cnt1),tx,squeeze(sdf(tinf{sel(jt)}(1),tinf{sel(jt)}(zt),:)),'r','LineWidth',3);
                        axis('tight');
                        axis('off');
                        %axis(ax5(cnt1),'off');
                        %axis(ax5(cnt1),'tight');axis(ax5(cnt1),'off');
                        [v,ix] = max(squeeze(sdf(tinf{sel(jt)}(1),tinf{sel(jt)}(zt),:)));
                        text(tx(ix),v+v/100*5,['Spikes/',num2str(pt(2)-pt(1)),'s: ',num2str(txt)],'parent',gca);
                        %text(tx(ix),v+v/100*5,['Spikes/',num2str(pt(2)-pt(1)),'s: ',num2str(txt)],'parent',ax5(cnt1));
                        xlim([-.2 1]);%xlim([-bt(1) pt(2)]);
                        %xlim(ax5(cnt1),[-bt(1) pt(2)]);
                    end;
                    
                    %%
                    Y =[];
                    Y(:,1) = squeeze(SC(tinf{sel(jt)}(1),:,1));
                    Y(:,2) = squeeze(SC(tinf{sel(jt)}(1),:,2));
                    plot(ax6(jt),[0 size(SC,2)+1],[thr2{sel(jt)} thr2{sel(jt)}],'k','LineWidth',3);
                    errorbar(ax6(jt),1:size(SC,2),Y(:,1),zeros(1,size(Y,1)),Y(:,2),'s','Color',[.75 .75 .75]);
                    errorbar(ax6(jt),tinf{sel(jt)}(2:end),Y(tinf{sel(jt)}(2:end),1)',Y(tinf{sel(jt)}(2:end),2)','s','Color',[.9 0 0],'MarkerFaceColor',[.9 0 0]);
                    axis(ax6(jt),'tight');
                    xlabel(ax6(jt),'Picture Nr.');
                    ylabel(ax6(jt),'Nr. of spikes');
                    title(ax6(jt),[chan,': cluster # ',num2str(sel(jt))]);                    
                    
                    plot(ax1,SCORE{tinf{sel(jt)}(1)}(:,1),SCORE{tinf{sel(jt)}(1)}(:,2),'.','Color',C(tinf{sel(jt)}(1),:));
                    
                    X = spike_data.waveformtime;
                    Y = squeeze(spike_data.waveform{tinf{sel(jt)}(1)}(1,:,:));
                    plot(ax2(jt),X,Y,'Color',C(tinf{sel(jt)}(1),:));
                    plot(ax2(jt),X,mean(Y,2),'Color',[1 1 1],'LineWidth',3);
                  
                    set(ax2(jt),'XTick',[min(spike_data.waveformtime) max(spike_data.waveformtime)]);
                    set(ax2(jt),'XTickLabel',round(get(ax2(jt),'XTick').*1e3*1e2)/1e2);
                    title(ax2(jt),['C#',num2str(sel(jt)),': ',num2str(length(spike_data.timestamp{tinf{sel(jt)}(1)}))]);
                    xlabel(ax2(jt),'Time [ms]');
                    ylabel(ax2(jt),'Amplitude [\muV]');
                    box(ax2(jt),'off');
                    axis(ax2(jt),'tight');
                    set(ax2(jt),'YTick',[min(get(gca,'YTick')) max(get(gca,'YTick'))]);
                    
                    bar(ax3(jt),dt_isi(dt_isi<=1500),(n_isi(tinf{sel(jt)}(1),dt_isi<=1500)./sum(n_isi(tinf{sel(jt)}(1),:))).*100,'FaceColor',C(tinf{sel(jt)}(1),:));
                    plot(ax3(jt),[min(dt_isi) max(dt_isi)],[1 1],'r--');
                    axis(ax3(jt),'tight');xlim(ax3(jt),[0 250]);set(ax3(jt),'XTick',[0 100 200]);
                    %set(gca,'XTickLabel',get(gca,'XTick').*1e3);
                    l = num2str(pct(tinf{sel(jt)}(1)));
                    title(ax3(jt),[l,' % < ',num2str(refP(1)*1e3),'ms']);
                    xlabel(ax3(jt),'ISI [ms]');
                    ylabel(ax3(jt),'Frequency [%]');
                    box(ax3(jt),'off');
                    
                    tv      = 10:-.05:1;%thres levels
                    
                    %create permutations
                    psel    = {};
                    psel2   = {};
                    for nt = 1:100
                        dum = randperm(size(SC,2));
                        dum(ismember(dum,tinf{sel(jt)}(2:end))) = [];
                        try
                            psel{nt} = dum(1:nrep);
                        catch
                             psel{nt} = dum;   
                        end;
                        psel2{nt} = setdiff(1:size(SC,2),[psel{nt} tinf{sel(jt)}(2:end)]);
                    end;
                    
                    
                    for ot = 2:length(tinf{sel(jt)})
                        
                        cnt2 = cnt2+1;
                        
                        TP  = zeros(1,length(tv));
                        FP  = zeros(1,length(tv));
                        TP2  = zeros(length(tv),100);
                        FP2  = zeros(length(tv),100);
                        
                        for zt = 1:length(tv)%loop over thresh
                            
                            if tv(zt)>1
                                ROCthr = mean(x1{tinf{sel(jt)}(1)}) + tv(zt)*std(x1{tinf{sel(jt)}(1)});
                            else
                                ROCthr = 0;
                            end;
                            
                            sc  = SC(tinf{sel(jt)}(1),:,1);
                            ix  = sc >= ROCthr;
                            sel2 = tinf{sel(jt)}(ot);
                            sel3 = setdiff(1:length(sc),sel2);
                            
                            TP(zt) = sum(ix(sel2))/length(sel2);%true positives
                            FP(zt) = sum(ix(sel3))/length(sel3);%false positives
                            
                            for nt = 1:length(psel)%load permutations
                                TP2(zt,nt) = sum(ix(psel{nt}))/length(psel{nt});%keep thresh levels constant
                                FP2(zt,nt) = sum(ix(psel2{nt}))/length(psel2{nt});
                            end;
                            
                        end;
                        
                        AUC1 = trapz(FP,TP);%AUC
                        
                        % calculate AUC for permuted images
                        AUC2 = zeros(1,size(FP2,2));
                        for nt = 1:size(FP2,2)
                            AUC2(nt) = trapz(FP2(:,nt),TP2(:,nt));
                        end;
                        p = (sum(AUC2 >= AUC1)/length(AUC2));
                        
                        for nt = 1:size(FP2,2)
                            plot(ax7(cnt2),FP2(:,nt),TP2(:,nt),'Color',[.75 .75 .75]);
                        end;
                        plot(ax7(cnt2),[0 1],[0 1],'r--');
                        plot(ax7(cnt2),FP,TP,'b*-');
                        title(ax7(cnt2),['AUC: ',num2str(AUC1),' ,p= ',num2str(p)]);
                        xlabel(ax7(cnt2),'False positive rate');
                        ylabel(ax7(cnt2),'True positive rate');
                    end;
                    
                end;
            end;
        end;
        tinf = {};
    end;
end;
return;
