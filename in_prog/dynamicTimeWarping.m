
%%
Y = mean(dum.wavf_decorr(find(dum.assignedCluster==1),:),1);%n1(1,:);
Yp = dum.wavf_decorr(2,:);

%%
D = zeros(length(Y),length(Yp))+Inf;
for it = 1:length(Y)
    
    for jt = 1:length(Yp)
        
        cost = sqrt((Y(it)-Yp(jt)).^2);
        D(it,jt) = cost;
    end;
    
end;

AC = zeros(size(D));
AC(1,2:end) = cumsum(D(1,2:end));
AC(2:end,1) = cumsum(D(2:end,1));

for it = 2:size(D,1)
    for jt = 2:size(D,2)
        AC(it,jt) = min([AC(it-1,jt-1) AC(it-1,jt) AC(it,jt-1)]) + D(it,jt);
    end;
end;

path = [];
d = [];
it = length(Y);
jt = length(Yp);
path = [jt it];
d = [jt it];
while it > 2 && jt > 2
    d = d-1;
    if it ==1
        jt = jt-1;
    elseif jt ==1
        it = it-1;
    else
        if AC(it-1,jt) == min([AC(it-1,jt-1) AC(it-1,jt) AC(it,jt-1)]) %&& max(abs([it-1 jt]-d)) < 50
            it = it-1;
        elseif AC(it,jt-1) == min([AC(it-1,jt-1) AC(it-1,jt) AC(it,jt-1)]) %&& max(abs([it jt-1]-d)) < 50
            jt = jt-1;
        else
            it = it-1;
            jt = jt-1;
        end;
    end;
    path = [ path;jt it];
end;
path = [path;1 1];

figure;
subplot(3,3,[1 4]);
hold on;
plot(Y,1:length(Y));
%plot(Y(path(:,1)),1:length(Y),'r');
axis tight;
set(gca,'XDir','reverse');
subplot(3,3,[8 9]);
hold on;
plot(1:length(Yp),Yp);
ix = flipud(path(:,2));
ix(find(diff(ix)==0))=[];
plot(1:length(Y),Yp(path(ix,2)),'r');
axis tight;
subplot(3,3,[2 3 5 6]);
imagesc((AC));
axis xy;
hold on;
plot([1 size(D,1)],[1 size(D,2)],'w');
for it = 2:size(path,1)    
    plot([path(it-1,1) path(it,1)],[path(it-1,2) path(it,2)],'gs-');    
end;

%%
ix1 = zeros(size(D,1),1);
ix2 = zeros(size(D,1),1);
ix1(1) = 1;
ix2(1) = 1;
ix1(end) = size(D,2);
ix2(end) = size(D,2);
ref = [2 2];
for it = 2:size(D,2)-1
    
    step = [];    
    step(1,:) = [0 0];
    step(2,:) = [1 0];
    step(3,:) = [0 1];
    step(4,:) = [1 1];
    
    [v,ix] = min(D(ref(1),ref(2))+[D(it,it) D(it+1,it) D(it,it+1) D(it+1,it+1)]);
    
    ix1(it) = it+step(ix,2);
    ix2(it) = it+step(ix,1);
    
    ref = [ix1(it) ix2(it)];%[ref(1)+step(ix,2) ref(2)+step(ix,1)];
    
end;

d = median(median(D));
figure;
subplot(3,3,[1 4]);
hold on;
plot(Y,1:length(Y),'r');axis tight;
plot(Y(ix1),1:length(Yp),'b');axis tight;
set(gca,'XDir','reverse');
subplot(3,3,[8 9]);
hold on;
plot(1:length(Yp),Yp,'b');axis tight;
plot(1:length(Yp),Yp(ix2),'r');axis tight;

subplot(3,3,[2 3 5 6]);
hold on;
imagesc(D);axis xy;%caxis([0 d])
axis tight;
for it = 2:length(ix1)
    plot([it-1 it],[it-1 it],'r-');
    plot(it,it,'ro');
    plot([ix1(it-1) ix1(it)],[ix2(it-1) ix2(it)],'g-');
    plot(ix1(it),ix2(it),'go');
end;

