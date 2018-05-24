function makePrinComp_plot(pcs,sel_ix,wvf,clu_ix)

ix = [];
for it = 1:length(sel_ix)
    if it ==1
        ix(it,:) = [1 sel_ix(it)];
    else
        ix(it,:) = [ix(it-1,2)+1 ix(it-1,2)+sel_ix(it)];
    end;
end;
sel_ix = ix;

c = .1:.1:1;
C = zeros(length(c)^3,3);
k = 0;
for it = 1:length(c)
    for jt = 1:length(c)
        for kt = 1:length(c)
            k = k+1;
            C(k,:) = [c(it) c(jt) c(kt)];
        end;
    end;
end;
[i1,~] = ind2sub(size(C),find(diff(abs(diff(C,[],2)),2)==0));
C(i1,:) = [];

%%
n = size(sel_ix,1);
ix = 1:n*3;
ix(3:3:end) = [];

ix2 = 1:n*3;
ix2 = setdiff(ix2,ix);

figure;
subplot(n,3,ix);
hold on;
d_ix = setdiff(1:size(sel_ix,1),clu_ix);
for it = 1:length(d_ix)
    plot3(pcs(sel_ix(d_ix(it),1):sel_ix(d_ix(it),2),1),pcs(sel_ix(d_ix(it),1):sel_ix(d_ix(it),2),2),pcs(sel_ix(d_ix(it),1):sel_ix(d_ix(it),2),3),'.','Color',[.75 .75 .75]);
end;

mem = [];
c = [];
for it =1:length(clu_ix)
    r_ix = randperm(size(C,1));
    r_ix = r_ix(1);
    
    if it > 1
    f =0;
    while f <1
    
        d = abs(diff([r_ix mem]));
        if d > 20
            f=1;
        end;
        r_ix = randperm(size(C,1));
        r_ix = r_ix(1);
    end;
    end;
    
    h = plot3(pcs(sel_ix(clu_ix(it),1):sel_ix(clu_ix(it),2),1),pcs(sel_ix(clu_ix(it),1):sel_ix(clu_ix(it),2),2),pcs(sel_ix(clu_ix(it),1):sel_ix(clu_ix(it),2),3),'.');
    set(h,'Color',C(r_ix,:));
    mem = r_ix;
    c(it,:) = C(r_ix,:);
end;

c2=0;
for it = 1:length(ix2)
    subplot(n,3,ix2(it));
    hold on;
    Y= squeeze(wvf{it})';
    if ismember(it,clu_ix)
        c2 = c2+1;
        plot(1:size(Y,2),Y,'Color',c(c2,:));
    else
        plot(1:size(Y,2),Y,'Color',[.75 .75 .75]);
    end;
    plot(1:size(Y,2),mean(Y,1),'k','LineWidth',3);
    axis tight;
    axis off;
end;