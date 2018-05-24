function image_resize(finf,iminf,sinf)
%%
if nargin == 0
    sinf.rsz = [360 720];
    sinf.fmt = 'bmp';
    iminf.ext = 'jpg';
    finf.p2d = '/home/rouxf/MATLAB/PTBcode/EMpairs_v4_2016-1007/image_data/Tune/';
    finf.fn = dir([finf.p2d,'*.',iminf.ext]);
end;
%%
if sinf.rsz(1) > sinf.rsz(2)
    sinf.rsz = fliplr(sinf.rsz);
end;
%%
for it = 1:length(finf.fn)
    
    fprintf([num2str(it),'/',num2str(length(finf.fn))]);
    
    [A] = imread([finf.p2d,finf.fn(it).name]);
    
    sz = size(A);
    
    if sum(sz(1:2) < sinf.rsz)==2
    else
        if sz(2) > sz(1)
            
            [B] = imresize(A,sinf.rsz);
            
        elseif sz(1) > sz(2)
            
            [B] = imresize(A, fliplr(sinf.rsz));
            
            
        end;
        sn = [finf.fn(it).name(1:end-4),'.',sinf.fmt];
        imwrite(B,[finf.p2d,sn],sinf.fmt);
        if ~strcmp(sn,finf.fn(it).name)
            delete([finf.p2d,finf.fn(it).name]);
        end;
    end;

    
    fprintf('\n');
end;