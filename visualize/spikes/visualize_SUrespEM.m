function visualize_SUrespEM(spike_dat,sel_ix1,sel_ix2,psth,cond,nLog,trl_all,dt)
%% make raster & spike density plot and bar plot
Fs = 1e3;

for it  = 1:length(sel_ix1)

    ts = spike_dat{sel_ix1(it)}.time{sel_ix2(it)+1}.*1e3;
    trl = spike_dat{sel_ix1(it)}.trial{sel_ix2(it)+1};
    trl = trl(find(ismember(trl,trl_all)));% retain trials during learning
    
    trlID = trl_all;
 
    figure;
    subplot(4,1,1:3);
    a =gca;
    hold on;
    plot([0 0],[0 nLog+1],'r');
    plot([2000 2000],[0 nLog+1],'r');
       
    x = [];
    for jt = 1:length(trlID)
        
        ix = find(trl == trlID(jt));        
        if isempty(ix)
            disp('toto');
        end;
        
        x = ts(ix);
        x(x<dt(1)) =[];
        x(x>dt(end)) =[];
        
        y = jt*ones(1,length(x));
        x = [x;x];
        y = [y-.5;y+.5];
        
%         ix1 = find(cond{it}(2,:)==1);
%         ix2 = find(cond{it}(2,:)==2);
        
%         if ismember(jt,cond{it}(1,ix1))
            line(x,y,'Color',[.5 .5 .5]);
%         elseif ismember(jt,cond{it}(1,ix2))
%             line(x,y,'Color','b');
%         else
%             line(x,y,'Color','r');
%         end;
        
    end;
    
    subplot(414);
    a = [a gca];
    Y = psth{it};
    
    M = sum(Y,1)/size(Y,1);
    SEM = std(Y,0,1)/sqrt(size(Y,1)-1);
    M = M/((dt(2)-dt(1))/Fs);
    SEM = SEM/((dt(2)-dt(1))/Fs);
    hold on;
    if any( SEM ~= 0 ) && any(~isnan( SEM ))
        plot([0 0],[min(M-SEM) max(M+SEM)],'r');
        plot([2000 2000],[min(M-SEM) max(M+SEM)],'r');
    else
        plot([0 0],[min(M) max(M)],'r');
        plot([2000 2000],[min(M) max(M)],'r');
    end;
    errorbar(dt,M,SEM,'s-','Color',[.25 .25 .25],'MarkerFaceColor',[.5 .5 .5]);
    axis tight;box on;
    %axis(a(1),'off');
    set(a,'XLim',[dt(1) dt(end)]);
    set(gcf,'Color','w');
    set(a(1),'YTick',[10:10:size(psth{it},1)]);
end;