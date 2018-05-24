function rasterPlot(ts,C,varargin)
if isempty(varargin)
    yl = 'Trial #';
else
    yl = varargin{1};
end;

hold on;
for it = 1:size(ts,1)    
    plot([ts{it,:}],it*ones(1,length([ts{it,:}])),'.','Color',C);
end;

ylabel(gca,yl);
xlabel(gca,'Time (s)');
axis(gca,'tight');

return;