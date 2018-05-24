function [sC] = getSphereCoords(s)

if nargin == 0
    %%
    s.x = 32;
    s.y = 41;
    s.z = 56;
    s.r = 5;
end;

gP = [s.x s.y s.z];

theta = linspace(0,2*pi);
phi = linspace(0,2*pi);
[theta,phi] = meshgrid(theta,phi);
[xs,ys,zs] = sph2cart(theta,phi,s.r);

sC.sX = (xs+gP(1));
sC.sY = (ys+gP(2));
sC.sZ = (zs+gP(3));
%%
function [sel] = sphere2gridPoint(sC,gC)

xb = [min(min(sC.sX)) max(max(sC.sX))];
yb = [min(min(sC.sY)) max(max(sC.sY))];
zb = [min(min(sC.sZ)) max(max(sC.sZ))];

k = 0;
sel = [];
for it = 1:size(gC,1)
    
    if sum(gC(it,:) >= [xb(1) yb(1) zb(1)])==3 && sum(gC(it,:) <= [xb(2) yb(2) zb(2)])==3
        k = k+1;
        sel(k) = it;
    end;
    
end;
