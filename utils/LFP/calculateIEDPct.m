function [pct,trlChck] = calculateIEDPct(ntrl,ix,pctTrsh,plt)
%%
trlChck = cell(1,size(ix,2));
pct = zeros(1,size(ix,2));

for it = 1:size(ix,2)
    
    x = ix(:,it);
    
    c = 0;    
    for jt = 1:length( x )
        if ~isempty( x{jt} )
            c = c+1;
            trlChck{it}(c) = jt;
        end;
    end;
    pct(it) = round(length([trlChck{it}])/ntrl*100)/100;
    
end;

if strcmp(plt,'y')
    figure;
    hold on;
    for it = 1:length( pct )
        if pct(it) <=pctTrsh
            plot(trlChck{it},it*ones(1,length(trlChck{it})),'ko','MarkerFaceColor','k');
        else
            plot(trlChck{it},it*ones(1,length(trlChck{it})),'ko','MarkerFaceColor','r');
        end;
        ylim(gca,[0 size(ix,2)+1]);
        set(gca,'YTick',1:size(ix,2));
        %set(gca,'YTickLabel',micDat.label);
        ylabel('Channel #');
        xlabel('Trial #');
    end;
end;