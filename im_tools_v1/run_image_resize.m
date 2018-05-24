%%
addpath('~rouxf/Desktop/im_tools_v1/');

sinf.rsz = [180 180];
sinf.fmt = 'jpg';
iminf.ext = 'jpg';
finf.p2d = '/media/rouxf/rds-share/Fred/cnEM_se06_2017-05-04/mix/';
finf.fn = dir([finf.p2d,'*.',iminf.ext]);
if strcmp(iminf.ext,'*')
    finf.fn(1:2) = [];
end;

image_resize(finf,iminf,sinf)
