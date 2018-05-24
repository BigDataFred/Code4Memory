function visualize_sorting(spdat)


wf = spdat.newSpikesNegative.*spdat.scalingFactor.*1e6;%convert to uV
%%
cid = spdat.useNegative;
n = zeros(length(cid),1);
ts_idx = cell(length(cid),1);
for it = 1:length(cid)
    n(it) = length(find(spdat.assignedNegative == cid(it)));
    ts_idx{it} = find(spdat.assignedNegative == cid(it));
end;

%%
[v,ix] = sort(n);
cid = cid(ix);
cid = flipud(cid);

sel_idx = ix(end);
del_idx = find(n< 100);

ts_idx(del_idx) = [];
cid(del_idx) = [];
n(del_idx)  = [];
%%
[pc,score,latent,tsquare] = princomp(wf);
%%
c = {};
c{1} = [0 0 0];
c{2} = [.75 .75 .75];
c{3} = [.9 0 0];%b
c{4} = [0 .9 0];%g
c{5} = [0 0 .9];%r
c{6} = [.9 .9 0];%y

c{7} = [.9 0 .9];%m
c{8} = [0 .9 .9];%c

 figure;

 h = zeros(length(cid),1);
 for it = 1:length(cid)
     
     sel = find(spdat.assignedNegative == cid(it));     
     f = 0;
     while f<1
         idx = randperm(length(c));
         idx = idx(1);
         if it >1
             if ismember(idx,h)
             else
                 h(it) = idx;
                 f =1;
             end;
         else
             h(it) = idx;
             f =1;
         end;
     end;
     subplot(222);
     a3 =gca;
     hold on;
     plot(score(sel,1),score(sel,2),'.','Color',c{idx});

     subplot(221);
     a1 =gca;
     hold on;
     plot(1:256,mean(wf(sel,:),1),'Color',c{idx});
%      subplot(223);
%      a2 =gca;
%      hold on;
%      ts = spdat.newTimestampsNegative(sel);
%      %ts = diff(ts)/1000;
%      plot(ts,it*ones(1,length(ts)),'.','Color',c{idx});
end;
axis(a1,'tight');
%axis(a2,'tight');
axis(a3,'tight');

%set(a2,'YLim',[0 length(cid)]);

allInds = [];
for it = 1:length(cid)
    allInds = [allInds find(spdat.assignedNegative == cid(it))];
end;

yl(1,:) = [0.9*min(score(allInds,1)) 1.1*max(score(allInds,1))];
yl(2,:) = [0.9*min(score(allInds,2)) 1.1*max(score(allInds,2))];
set(a3,'XLim',yl(1,:));
set(a3,'YLim',yl(2,:));
xlabel(a3,'PC1');
ylabel(a3,'PC2');

set(a1,'XTickLabel',round((get(a1,'XTick')/(32000*4).*1000)*100)/100);

ylabel(a1,'[\muV]');
xlabel(a1,'Time [ms]');

title(a1,'Average waveforms');
title(a3,'PCA');

set(gcf,'Color','w');

return;





