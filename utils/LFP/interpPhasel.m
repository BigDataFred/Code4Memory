function [phWBi,phNBi,ixTWB,ixTNB,ixPhWB,ixPhNB] = interpPhasel(ft,ft2)
%ft = wideband
%ft2 = narrow-band
phWBi   = [];
phNBi   = [];
ixTWB   = [];
ixTNB   = [];
ixPhWB  = [];
ixPhNB  = [];

if ~isempty(ft)
    [ixMin,ixMax] = findExtrma(ft);
end;

[ixMin2,ixMax2] = findExtrma(ft2);

if ~isempty(ft)
    selIx = zeros(length(ixMin2),2);
    for it = 1:length( ixMin2 )
        d = abs(ixMin-ixMin2(it));
        [v,ix] = min(d);
        selIx(it,1) = ix;
        d = abs(ixMax-ixMax2(it));
        [v,ix] = min(d);
        selIx(it,2) = ix;
    end;
    ixMin = ixMin(selIx(:,1));
    ixMax = ixMax(selIx(:,1));
    
    delIx = find(diff(ixMin)==0);
    ixMin(delIx) = [];
    ixMax(delIx) = [];
    
    ixTWB = ixMin;   
end;
ixTNB = ixMin2;

%% estimate ascending and descending phase for wide-band signal
if ~isempty(ft)
    [aIx,dIx] = interpPhaseQuadrants(ft,ixMin,ixMax);
    [delIx] = find(sum(isnan([aIx' dIx']),2)>0);
    
    ix = [ixMin aIx' ixMax dIx'];
    ix(delIx,:) = [];
    
    phWBi = repmat([0 90 180 270],[size(ix,1) 1]);
    for kt = 1:size(phWBi,1)
        phWBi(kt,:) = phWBi(kt,:)+(kt-1)*360;
    end;
    n = size(ix,1);
    
    ix = ix';       ix = ix(:);
    phWBi = phWBi';     phWBi = phWBi(:);
    
    phWBi = interp1(ix,phWBi,ix(1):ix(end));
    ixPhWB = ix(1):ix(end);
    
    ix = reshape(ix',[length(ix)/n n])';
    %phi = reshape(phi',[length(phi)/n n])';
    
    for kt = 1:n
        if kt <n
            selIx = find(ismember(ix(1):ix(end),ix(kt,1):ix(kt+1,1)-1));
        else
            selIx = find(ismember(ix(1):ix(end),ix(kt,1):ix(end)));
        end;
        phWBi(selIx) = phWBi(selIx)-(360*(kt-1));
    end;
    ix = ix';       ix = ix(:);
end;

%% estimate ascending and descending phase for narrow-band signal
[aIx,dIx] = interpPhaseQuadrants(ft2,ixMin2,ixMax2);
[delIx] = find(sum(isnan([aIx' dIx']),2)>0);

ix = [ixMin2 aIx' ixMax2 dIx'];
ix(delIx,:) = [];

phNBi = repmat([0 90 180 270],[size(ix,1) 1]);
for kt = 1:size(phNBi,1)
    phNBi(kt,:) = phNBi(kt,:)+(kt-1)*360;
end;
n = size(ix,1);

ix = ix';       ix = ix(:);
phNBi = phNBi';     phNBi = phNBi(:);

phNBi = interp1(ix,phNBi,ix(1):ix(end));
ixPhNB = ix(1):ix(end);

ix = reshape(ix',[length(ix)/n n])';
%phi = reshape(phi',[length(phi)/n n])';

for kt = 1:n
    if kt <n
        selIx = find(ismember(ix(1):ix(end),ix(kt,1):ix(kt+1,1)-1));
    else
        selIx = find(ismember(ix(1):ix(end),ix(kt,1):ix(end)));
    end;    
    phNBi(selIx) = phNBi(selIx)-(360*(kt-1));
end;
ix = ix';       ix = ix(:);

% figure;
% hold on;
% plot(x,'k','LineWidth',3);hold on;
% plot(ft,'r','LineWidth',3);
% %plot(ix1,s(ix1),'bo','MarkerFaceColor','b');
% %plot(ix2,s(ix2),'ro','MarkerFaceColor','r');
% plot(ix1*ones(1,2),[min(x) max(x)],'m--');
% plot(ix3*ones(1,2),[min(x) max(x)],'b--');
% axis tight;

% figure;
% subplot(311);
% hold on;
% plot(ix(1):ix(end),x(ix(1,1):ix(end,end)),'k');
% plot(ix(1):ix(end),ft(ix(1,1):ix(end,end)),'r');xlim([2501 3501]);set(gca,'XTick',sort([ix1' ix3']));
% set(gca,'XTickLabel',{});
% subplot(312);hold on;[ax,h1,h2] = plotyy(ix(1,1):ix(end,end),phi2,ix(1,1):ix(end,end),phi2+pi);
% set(ax,'Xlim',[2501 3501]);set(ax,'XTick',sort([ix1' ix3']));set(ax,'YTick',[-pi 0 pi]);
% ylim(ax(1),[-pi pi]);ylim(ax(2),[0 2*pi]);set(h2,'LineStyle','none','Color','r','Marker','o');
% set(ax(1),'YTick',[-pi 0 pi]);set(ax(2),'YTick',[0 pi pi*2]);
% set(ax(1),'YTickLabel',{'-\pi','0','\pi'});set(ax(2),'YTickLabel',{'0','\pi','2\pi'});
% set(ax,'XTickLabel',{});
% subplot(313);hold on;plot(ix(1,1):ix(end,end),phi3,'r');plot(ix(1,1):ix(end,end),yi,'k.');xlim([2501 3501]);
% set(gca,'XTick',sort([ix1' ix3']));set(gca,'YTick',[0 180 360]);
% set(gca,'XTickLabel',{});