%%
params3.trialave = 0;
step = size(ICs1{4},2)/length(trlSel);
ix = 1:step;

h1 = figure;
ax2 =[];
for kt = 1:size(ICs1{it},1)
    ax2(kt) = subplot(size(ICs1{it},1),1,kt);
    hold(ax2(kt),'on');
end;

h2 = figure;
ax1 = [];
for kt = 1:size(ICs1{it},1)
    ax1(kt) =subplot(size(ICs1{it},1),1,kt);
    hold(ax1(kt),'on');
end;
    
for jt = 1:length(trlSel)
            
    [Sf,fx] = mtspectrumc( gradient(ICs1{it}(:,ix))' ,params3);
    
    plot(ax2(1),fx,Sf(:,1),'r');
    plot(ax2(2),fx,Sf(:,2),'g');
    plot(ax2(3),fx,Sf(:,3),'b');
    plot(ax2(4),fx,Sf(:,4),'c');
    plot(ax2(5),fx,Sf(:,5),'m');
    plot(ax2(6),fx,Sf(:,6),'y');
    %plot(ax2(7),fx,Sf(:,7),'k');
    %plot(ax2(8),fx,Sf(:,8),'Color',[.5 .5 .5]);
    
    x = ICs1{it}(1,ix);
    ft = filtfilt(b1,1,[fliplr(x) x fliplr(x)]);
    ft = ft(length(ix)+1:length(ix)*2);
    plot(ax1(1),ICs1{it}(1,ix),'r');plot(ax1(1),ft,'k');axis tight;
    
    x = ICs1{it}(2,ix);
    ft = filtfilt(b1,1,[fliplr(x) x fliplr(x)]);
    ft = ft(length(ix)+1:length(ix)*2);
    plot(ax1(2),ICs1{it}(2,ix),'g');plot(ax1(2),ft,'k');axis tight;
    
    x = ICs1{it}(3,ix);
    ft = filtfilt(b1,1,[fliplr(x) x fliplr(x)]);
    ft = ft(length(ix)+1:length(ix)*2);
    plot(ax1(3),ICs1{it}(3,ix),'b');plot(ax1(3),ft,'k');axis tight;
    
    x = ICs1{it}(4,ix);
    ft = filtfilt(b1,1,[fliplr(x) x fliplr(x)]);
    ft = ft(length(ix)+1:length(ix)*2);
    plot(ax1(4),ICs1{it}(4,ix),'c');plot(ax1(4),ft,'k');axis tight;
    
    x = ICs1{it}(5,ix);
    ft = filtfilt(b1,1,[fliplr(x) x fliplr(x)]);
    ft = ft(length(ix)+1:length(ix)*2);
    plot(ax1(5),ICs1{it}(5,ix),'m');plot(ax1(5),ft,'k');axis tight;
    
    x = ICs1{it}(6,ix);
    ft = filtfilt(b1,1,[fliplr(x) x fliplr(x)]);
    ft = ft(length(ix)+1:length(ix)*2);
    plot(ax1(6),ICs1{it}(6,ix),'y');plot(ax1(6),ft,'k');axis tight;
%     
%     x = ICs1{it}(7,ix);
%     ft = filtfilt(b1,1,[fliplr(x) x fliplr(x)]);
%     ft = ft(length(ix)+1:length(ix)*2);
%     plot(ax1(7),ICs1{it}(7,ix),'k');plot(ax1(7),ft,'k');axis tight;
%     
%     x = ICs1{it}(8,ix);
%     ft = filtfilt(b1,1,[fliplr(x) x fliplr(x)]);
%     ft = ft(length(ix)+1:length(ix)*2);
%     plot(ax1(8),ICs1{it}(8,ix),'Color',[.5 .5 .5]);plot(ax1(8),ft,'k');axis tight;
    
    pause;
    for kt = 1:length(ax1);cla(ax1(kt));cla(ax2(kt));end;
    ix = ix+step;
    
end;