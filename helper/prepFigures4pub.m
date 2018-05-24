%%
nImages =10;
for it = 1:nImages
        
    set(it,'Color','w');
    set(it,'PaperPositionMode','auto');
    a = get(it,'Children');
    h = findobj(it,'Type','Line');
    for jt = 1:length(h)
        set(h(jt),'LineWidth',3);
    end;
    
    set(a,'Fontsize',14);
    set(a,'FontName','Arial');
    set(a,'LineWidth',3);
    %axis(a,'tight');
    
    for jt  = 1:length(a)
        set(get(a(jt),'XLabel'),'Fontsize',16)
        set(get(a(jt),'YLabel'),'Fontsize',16)
         set(a(jt),'TickLength',get(a(jt),'TickLength').*2);
         yl = get(a(jt),'YTick');
         %b = get(a(jt),'Xaxis');         
         %set(a(jt),'XTick',b.Limits);
         
         if min(yl)<0 && max(yl)>0
            set(a(jt),'YTick',[min(yl) 0 max(yl)]);
         else
            set(a(jt),'YTick',[min(yl) max(yl)]); 
         end;
    end;
    
end;
    
%%
savepath = '/home/rouxf/tmp/figuresIcon/';
set(gcf,'PaperPositionMode','auto');
for it = 10%1:nImages
    print(it,'-r800','-dtiff',[savepath,['pooledAnalysis',num2str(5),'.tif']]);
end;