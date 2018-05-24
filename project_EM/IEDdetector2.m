function [delIx,pct] = IEDdetector2(ERPimg)
delIx = [];
% 
% x = zeros(size(ERPimg,1),size(ERPimg,3));
% for it = 1:size(ERPimg,2)
%     dum = squeeze(ERPimg(:,it,:));
%     dum = dum./max(max(max(dum)));
%     x = x + dum;
% end;
% x = x./it;
% x = x./max(max(max(abs(ERPimg))));
% ix = find(abs(x)>=.85);
% 
% M = repmat(mean(x,2),[1 size(x,2)]);
% SD = repmat(std(x,0,2),[1 size(x,2)]);
% z = (x-M)./SD;
% 
% [i,j]  = ind2sub(size(x),ix);
% ix = ix(abs(z(ix))>4);
% 
% delIx = unique(i);
% if size(delIx,1) > size(delIx,2)
%     delIx = delIx';
% end;

delIx2 = {};
for it = 1:size(ERPimg,2)
    
    x =  squeeze(ERPimg(:,it,:));    
    M = repmat(mean(x,2),[1 size(x,2)]);
    SD = repmat(std(x,0,2),[1 size(x,2)]);
    z = (x-M)./SD;    
    
    x = x./max(max(x));    
    
    thr = abs(x)>=.85;
    
%     figure;
%     imagesc(thr);caxis([-1 1]);
    
    ix = find(thr == 1);
    %ix = ix(abs(z(ix))>median(iqr(abs(z')))*5);
    ix = ix(abs(z(ix))>6);
    
    [i,j]  = ind2sub(size(x),ix);
    y = unique(i)';
    delIx2{it} = y;
        
end;
delIx = delIx2;

pct = [];
for it = 1:length(delIx2)
    pct(it) = length([delIx2{it}])/size(ERPimg,1);
end;

% delIx2 = unique([delIx2{:}]);
% 
% if size(delIx2,1) > size(delIx2,2)
%     delIx2 = delIx2';
% end;
% 
% delIx = unique([delIx delIx2]);
return;

