function [ixMin,ixMax] = findExtrma(x)

% find local minima
ix1 = find(diff(sign(diff(x)))==2)+2;
ix1 = ix1(sign(x(ix1))==-1);

% find local maxima
ix3 = find(diff(sign(diff(x)))==-2)+2;
ix3 = ix3(sign(x(ix3))==1);

% prune peaks and troughs
dum = [];
dum(ix1) = -1;
dum(ix3) = 1;
chck = [find(dum~=0)' dum(dum~=0)'];
chck = [ chck [2;diff(chck(:,2),[],1)] ];

delIx = [];
p = find(chck(:,3)==0);
while ~isempty(p)
    ix = min(p);
    ix = [ix-1 ix];
    f = 0;
    while f<1
        if (ix(end)+1<=size(chck,1)) && (chck(ix(end)+1,3) == 0)
            ix = [ix ix(end)+1];
        else
            f =1;
        end;
    end;
    p(find(ismember(p,ix))) = [];
    
    [~,mIx] = max(abs(x(chck(ix,1))));
    del = setdiff(1:length(ix),mIx);
    delIx = [delIx;ix(del)'];
end;
chck(delIx,:) = [];
ix1 = chck(find(chck(:,2)==-1),1);
ix3 = chck(find(chck(:,2)==1),1);

if ix3(1)<ix1(1)
    ix3(1) = [];
end;
if ix1(end)>ix3(end)
    ix1(end) = [];
end;

ixMin = ix1;
ixMax = ix3;