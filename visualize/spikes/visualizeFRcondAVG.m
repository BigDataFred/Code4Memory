function [h] = visualizeFRcondAVG(sel_ix1,psth,cond,dt,bw,toi)


for it = 1:length(sel_ix1)
    
    T = dt(2)-dt(1);
    tInt = [toi(1):bw:toi(2)];
    x = [];FR = [];
    for ft = 1:length(tInt)-1
        nSpk = sum(psth{it}(:,find(dt >=tInt(ft) & dt < tInt(ft+1))),2);
        x(:,ft) = nSpk./T;% poststim spike count for each trial, note this is not firing rate
    end;
    FR = mean(x,2);
    
    
    condID = unique(cond{it}(2,:));
    
    figure;
    hold on;
    h = [];
    for jt = 1:length( condID )
                
        condIdx = find( cond{it}(2,:) == condID(jt));
        M = mean(FR(condIdx));
        SEM = std(FR(condIdx)/sqrt(size(FR,1)-1));
        QT = quantile(FR(condIdx),[.25 .75]);
        h(jt) = bar(jt,M);
        plot([jt jt],[M-SEM M+SEM],'r','LineWidth',3);
        %plot([jt-.25 jt+.25],[M M],'k','LineWidth',3);
        plot([jt-.25 jt+.25],[M+SEM M+SEM],'k','LineWidth',3);
        plot([jt-.25 jt+.25],[M-SEM M-SEM],'k','LineWidth',3);
        %plot([jt-.25 jt+.25],[QT(1) QT(1)],'b','LineWidth',3);
        %plot([jt-.25 jt+.25],[QT(2) QT(2)],'b','LineWidth',3);
        
        %plot(jt*ones(1,length(condIdx)),FR(condIdx),'ko','MarkerFaceColor','k');
        %plot([jt-.25 jt-.25],[QT(1) QT(2)],'b','LineWidth',3);
        %plot([jt+.25 jt+.25],[QT(1) QT(2)],'b','LineWidth',3);
    end;
    set(gca,'Fontsize',14);
    set(gcf,'Color','w');
    set(gca,'Xlim',[-1 length(condID)+2]);
    set(gca,'XTick',condID);
    set(gca,'XTickLabel','');
    
end;
