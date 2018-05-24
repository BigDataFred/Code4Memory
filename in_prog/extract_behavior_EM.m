function [RTs,acc,nBlocks] = extract_behavior_EM(rpath)
%%
if nargin == 0
    rpath = '/media/rouxf/rds-share/iEEG_DATA/MICRO/P02/fvSpEM/';
end;
%% check for existing sessions
chck = dir(rpath);
chck(1:2) = [];
ix = [];
c =0;
for it = 1:length(chck)
    if chck(it).isdir
        c =c+1;
        ix(c) = it;
    end;
end;

sesh = chck(ix);

if isempty(sesh)
    error('there is no corresponding data');
end;

%% loop over sessions and extract behavioral data
acc = zeros(length(sesh),3);
nBlocks = cell(2,length(sesh));
for nt = 1:length(sesh)
    
    p2d = [rpath,sesh(nt).name,filesep,'log_dat',filesep];
    fn = dir([p2d,'*LogDat.mat']);
    
    load([p2d,fn.name]);% load the data
    
    acc(nt,1) = length(ix{4}); % number of correct selections (100)
    acc(nt,2) = length(ix{5}); % number of 50/50 selections
    acc(nt,3) = length(ix{6}); % number of in-correct selections (100)
    
    ntrl = sum(acc(nt,:)); % number of trials
    
    acc(nt,1) = acc(nt,1)./ntrl; % pct correct
    acc(nt,2) = acc(nt,2)./ntrl; % pct 50/50
    acc(nt,3) = acc(nt,3)./ntrl; % pct incorrect
    
    nBlocks{1,nt} = length(diff(LogDat1.idx,[],2));
    nBlocks{2,nt} = diff(LogDat1.idx,[],2)';
    
    dum{1} = [RTs(ix{4})];
    dum{2} = [RTs(ix{5})];
    dum{3} = [RTs(ix{6})];
    
    RTs = dum;
    
end;

% %% calc stats
% m = mean(n,2);
% sem = std(n,0,2)/sqrt(size(n,2)-1);
% 
% h = zeros(1,length(m));
% figure;
% hold on;
% for it = 1:length(m)
%     h(it) = makeSEMbarGraph(it,m(it),sem(it));
% end;
% 
% set(h,'FaceColor',[.75 .75 .75]);
% set(h,'LineWidth',3);
% 
% %%
% function [h] = makeSEMbarGraph(x,m,sem)
% 
% h = bar(x,m);
% plot([x-.1 x+.1],[m+sem m+sem],'k','LineWidth',3);
% plot([x-.1 x+.1],[m-sem m-sem],'k','LineWidth',3);
% plot([x x],[m-sem m+sem],'k','LineWidth',3);


