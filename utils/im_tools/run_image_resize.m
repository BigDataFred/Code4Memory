%%
sinf.rsz = [360 720];
sinf.fmt = 'bmp';
iminf.ext = '*';
finf.p2d = '/home/adf/rouxf/Desktop/tuning_p04_201016_SE06/';
finf.fn = dir([finf.p2d,'*.',iminf.ext]);
if strcmp(iminf.ext,'*')
    finf.fn(1:2) = [];
end;

image_resize(finf,iminf,sinf)
