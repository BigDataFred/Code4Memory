%%
x = 1:100;
a = 1:10;

figure;
ax1=gca;
hold on;
title('orig spectrum');

figure;
ax2=gca;
hold on;
title('orig log-log');

figure;
ax3=gca;
hold on;
title('pred log-log');

figure;
ax4=gca;
hold on;
title('pred spectrum');

b  = [];
yp = [];
for it = 1:length(a)
    
    %subplot(211);
    fx = 1./(x.^a(it));
    plot(ax1,x,fx);
    
%     subplot(212);
     fx = 1./(x.^a(it));
     plot(ax2,log(x),log(fx));
    
     b(it,:) = regress(log(fx)',[ones(length(x),1) log(x)']);
     
     yp(:,it) = b(it,1)+b(it,2)*log(x)';
     
     plot(ax3,log(x),yp(:,it));
     
     plot(ax4,x,1./(x.^abs(b(it,2))));
end;
