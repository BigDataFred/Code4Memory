%%
addpath('/data/home1/froux/fieldtrip-20140304/');
ft_defaults;
path2files = '/data/home1/froux/DATASETS/DevelopmentData/figures/ICA_cleaning/';
path2files2 = '/data/projects2/Development/DevelopmentData/matFiles/WorkingMemory/ICA/';
%%
cfg = [];
cfg.viewmode = 'component';
cfg.layout = 'CTF275.lay';
%%
files1 = dir([path2files,'*_comp_idx1.fig']);
ekg_N = zeros(1,length(files1));
ekg_Ac = NaN(length(files1),20);
ekg_Fr = NaN(length(files1),1);
for it = 1:length(files1)
    
    cfile = dir([path2files2,files1(it).name(1:5),'_ICAcomps_*-Jun-2014.mat']);    
    load([path2files2,cfile.name]);
    
    ft_databrowser(cfg,comp);
    
    open([path2files,files1(it).name]);
    clc;
    fprintf(['processing :',files1(it).name(1:5),'\n']);
    k = 0;
    while k <1
        [s1] = input('Enter number of EKG components [n]','s');
        if ~isnan(str2double(s1))
            k = 1;
        end;
    end;
    ekg_N(it) = str2double(s1);
    
    for jt = 1:str2double(s1)
        [s2] = input(['Enter [0/1] for incorrect/correct detection of EKG component #', num2str(jt)],'s');
        ekg_Ac(it,jt) = str2double(s2);
    end;
    
    [s3] = input(['Enter number of falsely rejected EKG components #'],'s');
    ekg_Fr(it) = str2double(s3);
    
    close all;clc;
end;
% %%
% files2 = dir([path2files,'*_comp_idx2.fig']);
% eog_N = zeros(1,length(files2));
% eog_Ac = NaN(length(files2),4);
% for it = 1:length(files2)
%     
%     cfile = dir([path2files2,files2(it).name(1:5),'_DEVELOPMEN_ICAcomps_*-Apr-2014.mat']);
%     load([path2files2,cfile.name]);
%     
%     ft_databrowser(cfg,comp);
%     
%     open([path2files,files2(it).name]);
%     clc;
%     fprintf(['processing :',files2(it).name(1:5),'\n']);
%     
%     k = 0;
%     while k <1
%         [s1] = input('Enter number of EOG components [n]','s');
%         if ~isnan(str2double(s1))
%             k = 1;
%         end;        
%     end;
%     eog_N(it) = str2double(s1);
%     
%     for jt = 1:str2double(s1)
%         [s2] = input(['Enter [0/1] for incorrect/correct detection of EOG component #', num2str(jt)],'s');    
%         eog_Ac(it,jt) = str2double(s2);
%     end;
%     close all;clc;
%     
% end;
% eog_Ac(:,max(eog_N)+1:end) = [];
%%
files3 = dir([path2files,'*_comp_idx5.fig']);
emg_N = zeros(1,length(files3));
emg_Ac = NaN(length(files3),4);
for it = 1:length(files3)
    
    cfile = dir([path2files2,files3(it).name(1:5),'_DEVELOPMEN_ICAcomps_*-Apr-2014.mat']);
    load([path2files2,cfile.name]);
    
    ft_databrowser(cfg,comp);
    
    open([path2files,files3(it).name]);
    clc;
    fprintf(['processing :',files3(it).name(1:5),'\n']);
    k = 0;
    while k <1
        [s1] = input('Enter number of EMG components [n]','s');    
        if ~isnan(str2double(s1))
            k = k+1;
        end;
    end;
    emg_N(it) = str2double(s1);
    
    for jt = 1:str2double(s1)
        [s2] = input(['Enter [0/1] for incorrect/correct detection of EMG component #', num2str(jt)],'s');    
        emg_Ac(it,jt) = str2double(s2);
    end;
    close all;clc;
    
end;
emg_Ac(:,max(emg_N)+1:end) = [];
%%
files4 = dir([path2files,'*_comp_idx3_2.fig']);
files4b = dir([path2files,'*_comp_idx3_3.fig']);
ssp_N = zeros(1,length(files4));
ssp_Ac = zeros(1,length(files4));
for it = 1:length(files4)
    
    open([path2files,files4(it).name]);
    open([path2files,files4b(it).name]);
    clc;
    fprintf(['processing :',files4(it).name(1:5),'\n']);
    
    [s1] = input('Enter number of SSP components [n]','s');    
    ssp_N(it) = str2double(s1);
    
    [s2] = input('Enter [0/1] for incorrect/correct detection','s');    
    ssp_Ac(it) = str2double(s2);
    
    close all;clc;
    
end;
%%
files5 = dir([path2files,'*_comp_idx4_1.fig']);
files5b = dir([path2files,'*_comp_idx4_2.fig']);
nsc_N = zeros(1,length(files5));
nsc_Ac = zeros(1,length(files5));
for it = 1:length(files5)
    
    open([path2files,files5(it).name]);
    open([path2files,files5b(it).name]);
    clc;
    fprintf(['processing :',files5(it).name(1:5),'\n']);
    
    [s1] = input('Enter number of non-stationary components [n]','s');    
    nsc_N(it) = str2double(s1);
    
    [s2] = input('Enter [0/1] for incorrect/correct detection','s');    
    nsc_Ac(it) = str2double(s2);
    
    close all;clc;
    
end;
%%
path2files2 = '/data/projects2/Development/DevelopmentData/matFiles/WorkingMemory/rawSig/';
files6 = dir([path2files2,'*_DEVELOPMEN_basic_preproc_pre_AND_postICA_powerSpectra_29-May-2014.mat']);
SNR1 = zeros(length(files6),256);
SNR2 = zeros(length(files6),256);
SNR3 = zeros(length(files6),256);
for jt = 1:length(files6)
    fprintf([num2str(jt),'/',num2str(length(files6))]);
    load([path2files2,files6(jt).name]);
    
    SNR1(jt,:) = squeeze(mean(20.*log10(freq3.powspctrm./freq1.powspctrm),2));
    SNR2(jt,:) = squeeze(mean(20.*log10(freq3.powspctrm(:,freq3.freq<=20)./freq1.powspctrm(:,freq3.freq<=20)),2));
    SNR3(jt,:) = squeeze(mean(20.*log10(freq3.powspctrm(:,freq3.freq>=20)./freq1.powspctrm(:,freq3.freq>=20)),2));
    fprintf('\n');
end;
ix = find(squeeze(mean(SNR1,2))>10);
SNR1(ix,:,:) = [];
SNR2(ix,:,:) = [];
SNR3(ix,:,:) = [];
%%
h = [];
figure;
subplot(2,3,1);
hold on;
h(1) = plot(length(find(ekg_Ac==0))/length(files1),length(find(ekg_Ac==1))/length(files1),'ro','LineWidth',2);
h(2) = plot(length(find(eog_Ac==0))/length(files2),length(find(eog_Ac==1))/length(files2),'bo','LineWidth',2);
h(3) = plot(length(find(emg_Ac==0))/length(files3),length(find(emg_Ac==1))/length(files3),'ko','LineWidth',2);
h(4) = plot(length(find(ssp_Ac==0))/length(files4),length(find(ssp_Ac==1))/length(files4),'go','LineWidth',2);


legend(h,'EKG','EOG','EMG','SSP');

xlim([0 1]);
ylim([0 1]);
box on;grid on;
set(gca,'YTick',[0 .5 1]);
set(gca,'XTick',[0 .5 1]);
set(gca,'Fontsize',14);
xlabel('False detection [%]','Fontsize',14);
ylabel('Correct detection [%]','Fontsize',14);

subplot(2,3,2);
[n1,x1] = hist(ekg_N,[1:1:15]);
[n2,x2] = hist(eog_N,[1:1:15]);
[n3,x3] = hist(emg_N,[1:1:15]);
[n4,x4] = hist(ssp_N,[1:1:15]);

Nn1 = abs(sum(n1)-97);
Nn2 = abs(sum(n2)-97);
Nn3 = abs(sum(n3)-97);
Nn4 = abs(sum(n4)-97);

x1 = [0 x1];
x2 = [0 x2];
x3 = [0 x3];
x4 = [0 x4];

n1 = [Nn1 n1];
n2 = [Nn2 n2];
n3 = [Nn3 n3];
n4 = [Nn4 n4];

hold on;

Y = [n1' n2' n3' n4'];

c = {'b','r','k','g'};
for it = 1:size(Y,1)
    
    [y,ix] = sort(Y(it,:),'descend');
    for jt = 1:length(y)
        bar(it-1,y(jt),.5,c{ix(jt)});
    end;
end;

set(gca,'XTick',[0:2:15]);

xlim([0 max(x1)+1]);
ylim([0 100]);
set(gca,'YTick',[0 50 100]);
set(gca,'Fontsize',14);
xlabel('Number of components','Fontsize',14);
ylabel('Count','Fontsize',14);

subplot(233);
hold on;
bar([1],[mean(mean(SNR1))]);
bar([2],[mean(mean(SNR2))]);
bar([3],[mean(mean(SNR3))]);

plot([1 1],[mean(mean(SNR1)) mean(mean(SNR1)) - (std(mean(SNR1,2),0,1)./sqrt(size(SNR1,1)-1))],'k','LineWidth',3);
plot([2 2],[mean(mean(SNR2)) mean(mean(SNR2)) - (std(mean(SNR2,2),0,1)./sqrt(size(SNR2,1)-1))],'k','LineWidth',3);
plot([3 3],[mean(mean(SNR3)) mean(mean(SNR3)) - (std(mean(SNR3,2),0,1)./sqrt(size(SNR3,1)-1))],'k','LineWidth',3);

plot([.9 1.1],[mean(mean(SNR1)) - (std(mean(SNR1,2),0,1)./sqrt(size(SNR1,1)-1)) mean(mean(SNR1)) - (std(mean(SNR1,2),0,1)./sqrt(size(SNR1,1)-1))],'k','LineWidth',3);
plot([1.9 2.1],[mean(mean(SNR2)) - (std(mean(SNR2,2),0,1)./sqrt(size(SNR2,1)-1)) mean(mean(SNR2)) - (std(mean(SNR2,2),0,1)./sqrt(size(SNR2,1)-1))],'k','LineWidth',3);
plot([2.9 3.1],[mean(mean(SNR3)) - (std(mean(SNR3,2),0,1)./sqrt(size(SNR3,1)-1)) mean(mean(SNR3)) - (std(mean(SNR3,2),0,1)./sqrt(size(SNR3,1)-1))],'k','LineWidth',3);


h = findobj(gca,'Type','patch');
set(h(1),'FaceColor',[1 1 1]);
set(h(2),'FaceColor',[.75 .75 .75]);
set(h(3),'FaceColor',[0 0 0]);
xlim([0 4]);
set(gca,'XTick',[1:3]);
set(gca,'XTickLabel',{'all' '0-20Hz' '20-100Hz'},'Fontsize',14);
ylabel('SNR [dB]','Fontsize',14);
set(gca,'Fontsize',14);

cfg = [];
cfg.layout = 'CTF275.lay';
cfg.parameter = 'avg';
cfg.electrodes = 'off';
cfg.comment = 'no';
cfg.zlim = [-2 1];

subplot(234);
dum = struct;
dum.label = freq3.label;
dum.avg = squeeze(mean(SNR1,1))';
dum.time = 1;
dum.dimord = 'chan_time';
ft_topoplotER(cfg,dum);
cb = colorbar;
zlab = get(cb,'YLabel');
set(zlab,'String','SNR [dB]','Fontsize',14);
title('0-100Hz','Fontsize',14);
set(gca,'Fontsize',14);
subplot(235);
dum = struct;
dum.label = freq3.label;
dum.avg = squeeze(mean(SNR2,1))';
dum.time = 1;
dum.dimord = 'chan_time';
ft_topoplotER(cfg,dum);
cb = colorbar;
zlab = get(cb,'YLabel');
set(zlab,'String','SNR [dB]','Fontsize',14);
title('0-20Hz','Fontsize',14);
set(gca,'Fontsize',14);
subplot(236);
dum = struct;
dum.label = freq3.label;
dum.avg = squeeze(mean(SNR3,1))';
dum.time = 1;
dum.dimord = 'chan_time';

ft_topoplotER(cfg,dum);
cb = colorbar;
zlab = get(cb,'YLabel');
set(zlab,'String','SNR [dB]','Fontsize',14);
title('20-100Hz','Fontsize',14);
set(gca,'Fontsize',14);
set(gcf,'Color','w');
set(gcf,'PaperPositionMode','auto');
print(gcf,'-dtiff','-zbuffer',[pwd,'/ICAcleaning_summary_',date,'.tif'])