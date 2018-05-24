%%
imdir = '/home/adf/rouxf/Desktop/tuning_p05_se07_07122016/';
imf = dir([imdir,'*.jpg']);
%imf(1:2) = [];

savedir = [imdir,filesep,'resized',filesep];

chck = dir([imdir,filesep,'resized']);
if isempty(chck)
    mkdir([imdir,filesep,'resized']);
end;

%%
for it = 1:length(imf)
    
    fprintf([num2str(it),'/',num2str(length(imf))]);
    
    Xi = []; Yi = []; X = []; Y = [];
    
    [V,MAP,ALPHA] = imread([imdir,imf(it).name]);
    
    dim = size(V);
    [X,Y] = meshgrid(1:dim(2),1:dim(1));
    
    pct1 = ceil(dim(1)/1920*100);
    pct2 = ceil(dim(2)/1080*100);
    
    if pct1 <100
        n1 = (dim(1)+(((1920/100)*25)-dim(1)));
    else
        n1 = (1920/100*25);
    end;
    
    if pct2 <100
        n2 = (dim(2)+(((1920/100)*25)-dim(2)));
    else
        n2 = (1920/100*25);
    end;
    
    [Xi,Yi] = meshgrid(linspace(1,dim(2),n2),linspace(1,dim(1),n1));
    %%
    Yi1 = interp2(X,Y,double(squeeze(V(:,:,1))),Xi,Yi);
    Yi2 = interp2(X,Y,double(squeeze(V(:,:,2))),Xi,Yi);
    Yi3 = interp2(X,Y,double(squeeze(V(:,:,3))),Xi,Yi);
    
    %     PSF = fspecial('motion', 15);
    %     noise_mean = 0;
    %     noise_var = 0.0001;
    %
    %     estimated_nsr = noise_var / var(Yi1(:));
    %     Yi1 = edgetaper(Yi1, PSF);
    %     Yi1 = deconvwnr(Yi1, PSF, estimated_nsr);
    %
    %     estimated_nsr = noise_var / var(Yi2(:));
    %     Yi2 = edgetaper(Yi2, PSF);
    %     Yi2 = deconvwnr(Yi2, PSF, estimated_nsr);
    %
    %     estimated_nsr = noise_var / var(Yi3(:));
    %     Yi3 = edgetaper(Yi3, PSF);
    %     Yi3 = deconvwnr(Yi3, PSF, estimated_nsr);
    
    Yi1 = uint8(Yi1);
    Yi2 = uint8(Yi2);
    Yi3 = uint8(Yi3);
    
    Yi = cat(3,Yi1,Yi2,Yi3);
    
    [~,fn,ext] = fileparts(imf(it).name);
    fn(regexp(fn,'-')) = [];
    fn(regexp(fn,' ')) = [];
    fn(regexp(fn,'_')) = [];
    
    imwrite(Yi,[savedir,[fn,'_2'],ext],ext(2:end));
    
    fprintf('\n');
end;