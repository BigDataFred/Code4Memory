function visualizeSpkRaster(tx,spkRaster,pval,tx2,iFR,trsh,fr)
figure;
subplot(4,1,1:2);
hold on;
for it = 1:size(spkRaster,1)
    x = tx(spkRaster(it,:)~=0);
    y = it*ones(1,length(x));
    x = [x;x];
    y = [y-.5;y+.5];
    line(x,y,'Color','k');
end;
plot([0 0],[0 size(spkRaster,1)+1],'r');
plot([2e3 2e3],[0 size(spkRaster,1)+1],'r');
axis tight;
if ~isempty(pval)
    title(min(pval));
end;
subplot(4,1,3:4);
hold on;
plot(tx,iFR);
plot(tx2,fr);
plot([tx(1) tx(end)],[trsh trsh],'r--');
xlim([tx(1) tx(end)]);
%figure;plot(tx,iFR>trsh);ylim([-.5 2]);