function [emg_idx] = emg_comp_detect3(comp)
%%
cfg = [];
if length(comp.topolabel{end})==5
    cfg.layout = 'CTF275.lay';
else
    cfg.layout = 'neuromag306planar.lay';
end;

lay = ft_prepare_layout(cfg);
%%
chck = zeros(length(lay.label),1);
for it = 1:length(lay.label)
    chck(it) = any(strcmp(lay.label(it),comp.topolabel(:)));
end;

sel = find(chck==1);
%%
x = lay.pos(sel,1);
y = lay.pos(sel,2);
%%
xi = linspace(min(x),max(x),length(sel));
yi = linspace(min(y),max(y),length(sel));
%%
ix1 = [find(sign(xi-(-.2))==-1) find(sign(xi-.2)==1)];
ix2 = [find(sign(yi-(-.2))==-1) find(sign(yi-.2)==1)];
%%
D = zeros(length(comp.label),length(comp.topolabel),length(comp.topolabel));
for it = 1:length(comp.label)
    [Xi,Yi,Zi] = griddata(x',y,comp.topo(:,it),xi',yi,'v4');
    
    [Fx,Fy] = gradient(Zi);
    
    Zi = (abs(Fx)+abs(Fy));
    
    D2 = zeros(size(Zi));
    D2(ix1,:) = D2(ix1,:) + Zi(ix1,:);
    D2(:,ix2) = D2(:,ix2) + Zi(:,ix2);
        
    D(it,:,:) = D2;
end;
%%
M = repmat(mean(D,1),[size(D,1) 1 1]);
SD = repmat(std(D,[],1),[size(D,1) 1 1]);

Z = (D-M)./SD;
Z(isnan(Z))= 0;
%%
for it = 1:size(Z,1)
    d(it) = sum(sum(Z(it,:,:)));
end;
d = d-mean(d);
d = d./std(d);

[d,s_idx] = sort(d);
%%
emg_idx = s_idx(find(sign(d-0.01)==1));

