%%
clear;clc;

%%
restoredefaultpath;
addpath('~rouxf/fieldtrip-20161009/');
addpath(genpath('~rouxf/AnalysisFred/EM/'));
ft_defaults;

%%
p2d = '/media/rouxf/rds-share/iEEG_DATA/MICRO/P04/fvsPEM/';
d = dir([p2d,'2016-*']);

%%

dum = struct;
dum.trial = cell(1,24);
dum.time = cell(1,24);
dum.unit = cell(1,24);
dum.waveform = cell(1,24);
dum.timestamp = cell(1,24);

%% concatenate across runs
for jt = 1:24
    dat = {};
    for it =1:length(d)
        
        f = dir([p2d,d(it).name , filesep , 'spike_dat', filesep , '*_EMtask_sortedSPIKEdata.mat']);
        load([ p2d , d(it).name , filesep , 'spike_dat', filesep , f.name ]);
        
        dat = spike_data{2};
        
        ch = {};
        for zt = 1:length(dat.label)
            ch{zt} = dat.label{zt}(regexp(dat.label{zt},'_')+1:end);
        end;
        
        chck1 = regexp(ch,'L');
        chck2 = regexp(ch,'R');
        
        c = 0;sel_idx1 = [];
        for zt = 1:length(chck1)
            if ~isempty(chck1{zt})
                c = c+1;
                sel_idx1(c) = zt;
            end;
        end;
        
        c = 0;sel_idx2 = [];
        for zt = 1:length(chck2)
            if ~isempty(chck2{zt})
                c = c+1;
                sel_idx2(c) = zt;
            end;
        end;
        
        ch2 = [];
        for zt = 1:length(ch)
            ch2(zt) = str2double(ch{zt}(2:end));
        end;
        
        [~,s_idx1] = sort(ch2(sel_idx1));
        [~,s_idx2] = sort(ch2(sel_idx2));
        
        sel_idx1 = sel_idx1(s_idx1);
        sel_idx2 = sel_idx2(s_idx2);
        
        s_idx = [sel_idx1 sel_idx2]';
        
        
        dat.label = dat.label(s_idx);
        dat.trial = dat.trial(s_idx);
        dat.time = dat.time(s_idx);                
        dat.unit = dat.unit(s_idx);
        dat.waveform(s_idx);
        dat.timestamp(s_idx);
        
        dum.trial{jt} = [dum.trial{jt} ((it-1)*49)+dat.trial{jt}];
        dum.time{jt} = [dum.time{jt} dat.time{jt}];
        dum.unit{jt} = [dum.unit{jt} dat.unit{jt}];
        dum.waveform{jt} = [squeeze(dum.waveform{jt});squeeze(dat.waveform{jt})'];
        dum.timestamp{jt} = [dum.timestamp{jt} double(dat.timestamp{jt}-dat.hdr.FirstTimeStamp)/1e6];
    end;
    
end;
%%
r = 0.1:.1:1;
g = 0.1:.1:1;
b = 0.1:.1:1;

C = zeros(length(r)*length(g)*length(b),3);
c = 0;
for it = 1:length(r)
    for jt = 1:length(g)
        for zt = 1:length(b)
            c = c+1;
            C(c,:) = [r(it) g(jt) b(zt)];
        end;
    end;
end;

%%
addpath(genpath('~rouxf/chronux_2_11/'));

dt = [-.5:1e-3:6];
trl = unique([dum.trial{:}]);
figure;
draw = [];
for jt = 1%1:24
    %subplot(3,8,jt);
    %hold on;
    cID = unique(dum.unit{jt});
    
    S1 = []; R1 = []; Serr1 = [];
    figure;
    for zt = 8%1:length(cID)
                figure;
                hold on;
        
        f = 0;
        while f <1
            sel = randperm(size(C,1));
            sel = sel(1);
            if ~ismember(sel(1),draw)
                f = 1;
            end;
        end;
        draw(jt) = sel;
        
        spike_times1 = zeros(length(dt),length(trl));
        ix = 0; wvf = []; isi = [];
        for it = 1:length(trl)
            
            idx = find(dum.trial{jt} == trl(it) & dum.unit{jt} == cID(zt));
            
            if ~isempty(idx)
                ix = ix(end) + 1: ix(end) + length(idx);
                wvf(ix,:) = dum.waveform{jt}(idx,:);
                isi = [isi diff(dum.timestamp{jt}(dum.unit{jt} == cID(zt)))];
            end;
            
            x = dum.time{jt}(idx);
            x = x(x>=-.75 & x<=4);
            
            c = zeros(1,length(dt));
            if ~isempty(x)
                
                for gt = 1:length(dt)-1
                    
                    b = find(x >= dt(gt) & x < dt(gt+1));
                    c(gt) = length(b);
                end;
                                                                            
                y = trl(it)*ones(1,length(x));
%                 
                x = [x;x];
                y = [y-.5;y+.5];
                
                line(x,y,'Color','k');
            end;
            
            spike_times1(:,it) = c;
                                    
        end;
        
        params              = [];
        params.pad          = 0;
        params.fpass        = [2 25];
        params.tapers        = [1 1];
        params.err          = [1 0.001];
        params.trialave     = 1;
        params.Fs           = 1000;
        
        [S1(:,zt),f1,R1(:,zt),Serr1(:,:,zt)] = mtspectrumpb(spike_times1,params);
        
        plot([0 0],[min(trl) max(trl)],'r');
        plot([2 2],[min(trl) max(trl)],'r');
        
        axis tight;
        set(gca,'YTick',[min(trl) max(trl)]);
        set(gca,'XTick',[-.75 4]);
        xlim([-.75 4]);
        
        %subplot(ceil(length(cID)/5),5,zt);
        figure;
        plot(dat.waveformtime,wvf,'Color',C(sel,:));hold on;plot(dat.waveformtime,mean(wvf,1),'k','LineWidth',3);axis tight
        title(size(wvf,1));
        set(gca,'XTickLabel',get(gca,'XTick').*1e3);
        
        figure;
    hold on;
    plot(f1,squeeze(Serr1(:,:,zt))./(ones(2,size(S1,1)).*R1(:,zt)),'b');
    plot(f1,S1(:,zt)./(ones(size(S1(:,zt),1),1)*R1(:,zt)),'c');
    axis tight;
    end;
%     subplot(3,8,jt);
    
end;
%%
figure;
a = [];
for it = 1:length(psth)
    subplot(3,8,it);
    a(it) = gca;
    plot(psth.time,psth.avg);
    
end;
axis(a,'tight');