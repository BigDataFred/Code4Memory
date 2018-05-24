%% sanity check

Fs = 1000;
t = 0:1/Fs:6;

x = zeros(1,length(t));

x = x + sin(2*pi*14.*t);
x = x>=.999;

spike_times1 = x;

t = t.*1e3;
figure;
subplot(211);
plot(t,x,'b.');xlim([0 1e3]);

ts = t(find(x~=0));

spike_times3(1).times = ts';

dt = t(2)-t(1);
bins = min(t):dt:max(t);
c = zeros(length(bins),1);
for pt = 1:length(bins)-1
    
    if pt < length(bins)-1
        b = find(ts >= bins(pt) & ts < bins(pt+1));
    else
        b = find(ts >= bins(pt) & ts <= bins(pt+1));
    end;
    if ~isempty(b)
        c(pt) = c(pt) + 1;
    end;
    
end;

subplot(212);
plot(bins,c,'b.');xlim([0 1e3]);

spike_times2 = c;

figure;
plot(1:length(spike_times3(1).times),spike_times3(1).times,'b.');
ylim([0 1000]);

params              = [];
params.Fs           = 1000;
params.pad          = 0;
params.fpass        = [2 25];
params.tapers        = [3 5];
params.err          = [1 0.001];

%[S1,f1,R1,Serr1]=mtspectrumpb(spike_times1,params,0);
[S2,f2,R2,Serr2]=mtspectrumpb(spike_times2,params,0);

%[S3,f3,R3,Serr3]=mtspectrumpt(spike_times3,params);

figure;
hold on;
plot(f2,S2./R2,'r');
plot(f2,Serr2(1,:)./R2,'g-');
plot(f2,Serr2(2,:)./R2,'g-');

%plot(f1,S1./R1,'b')
%plot(f3,S3./R3,'k')
%%
clear;
clc;
addpath(genpath('~rouxf/chronux_2_11/'));

p2d = '/media/rouxf/rds-share/iEEG_DATA/MICRO/P04/fvsPEM';
d = dir([p2d,filesep,'2016-*']);
spike_all = {};%cell(1,length(d));
for it = 1:length(d)
    
    f = dir([p2d,filesep,d(it).name,filesep,'spike_dat',filesep,'04_fVSp__2016-*_LogFile_EMtask_sortedSPIKEdata.mat']);
    
    load([p2d,filesep,d(it).name,filesep,'spike_dat',filesep,f.name]);
    
    spike_all{it} = spike_data{2};
end;
%%
dum = struct;
for it = 1:length(spike_all{1}.label)
    dum.trial{it} = [];
    dum.time{it} = [];
    dum.unit{it} = [];
    for jt = 1:length(spike_all)
        dum.trial{it} = [dum.trial{it} (jt-1)*49+spike_all{jt}.trial{it}];
        dum.time{it} = [dum.time{it} spike_all{jt}.time{it}];
        dum.unit{it} = [dum.unit{it} spike_all{jt}.unit{it}];
    end;
end;
%%


a = [];
lm = [];
 figure;  
for zt = 1:24
    ts = dum.time{zt}';
    %ts = ts - spike_data{2}.hdr.FirstTimeStamp;
    %ts = double( ts./1e6 );
    %ts = double(ts).*1000;% convert s to ms
    
           c = dum.unit{zt};
    cID = unique(c);
    
    for vt = 1:length(cID)
    
    trl =  dum.trial{zt}';
    tID = unique(trl);%these are the trial indexes for each spike time
    
    %%
    spike_times1 = []; spike_times2 = []; sel_idx = [];n = []; k =0;
    
 
        for it = 1:length(tID)
            
            ix = [];
            ix = intersect(find(trl == tID(it)),find(c == cID(vt))); % find those spike times that correspond to a trial
            
            n(it) = length(ix);% number of spike times per trial
            
            %if ~isempty(ix) %length(find(ts(ix)>=2000 & ts(ix) <= 6000)~=0) >2% check if spike times are within time range
            k = k+1;%increment counter
            %[y,~] = hist(ts(ix),[2000:1:6000]);
            %[y,~] = hist(ts(ix),[2:0.001:6]);
            bins = 2:1e-3:6;
            c = zeros(length(bins),1);
            for pt = 1:length(bins)-1
                
                if pt < length(bins)-1
                    b = find(ts(ix) >= bins(pt) & ts(ix) < bins(pt+1));
                else
                    b = find(ts(ix) >= bins(pt) & ts(ix) <= bins(pt+1));
                end;
                if ~isempty(b)
                    c(pt) = c(pt) + length(b);
                end;
                
            end;
            spike_times1(:,k) = c;
            spike_times2(k).times = ts(ix);
            idx = find( ts(ix) >= 2 & ts(ix) <= 6);
            spike_times2(k).times = spike_times2(k).times(idx);
            sel_idx(k) = it;
            %end;
        end;
        %%
        zt
        params              = [];
        params.pad          = 0;
        params.fpass        = [2 25];
        params.tapers        = [1 1];
        params.err          = [1 0.001];
        params.trialave     = 1;
        subplot(3,8,zt);
        hold on;
        
        if ~isempty(spike_times1)
            params.Fs           = 1000;
            [S1,f1,R1,Serr1]=mtspectrumpb(spike_times1,params);
            %[S1,f1,R1,Serr1]=mtspectrumpb(spike_times1,params,0);
            Y1 = S1/R1;%mean(S1,2);%./(ones(size(S1,1),1)*R1'),2);%mean(S1,2)./mean(R1);%
            Err = Serr1./R1;%squeeze(mean(Serr1,3));%squeeze(mean(Serr1./shiftdim(repmat(R1,[1 size(Serr1,1) size(Serr1,2)]),1),3));
            %Serr1 = mean(Serr1./(ones(size(S1,1),1)*R1'),3);
            plot(f1,Y1,'COlor',[.9 0 0]);
            %plot(f1,Err(1,:),'r');
            %plot(f1,Err(2,:),'r');
            %title(num2str(mean(R1)));
        else
            disp('warning')
        end;
        % if ~isempty(spike_times2)
        %     params.Fs           = 32000;
        %     [S2,f2,R2,Serr2]=mtspectrumpt(spike_times2,params);
        %     %[S2,f2,R2,Serr2]=mtspectrumpt(spike_times2,params,0,[2000:1000:6000]);
        %     Y2 = mean(S2./(ones(size(S2,1),1)*R2),2);
        %     plot(f2,Y2,'r');
        %
        % end;
        axis tight;
        
        lm(zt,:) = [min(Y1) max(Y1)];% min(Y2) max(Y2)];
        a(zt) = gca;
    end;
end;
% for it = 1:length(a)
% ylim(a(it),[min(min(lm)) max(max(lm))]);
% end;