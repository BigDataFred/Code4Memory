%%
Fs = 1e3;
t = -1:1/Fs:3;

ix = find(t>0.5 & t<1);

[spk] = zeros(1,length( t ) );
spk(ix(1):25:ix(end)) = 1;

[xc,lag] = xcorr( spk-mean(spk) , 100, 'coeff' );

figure;
plot(t,spk);
axis tight;
ylim([-1 2]);

figure;
plot( lag,xc );
axis tight;

%%
xc = [];
ix = 1:100;
c=0;
for it = 1:length( spk )
    
    if (it>100) && (it < length( spk )-100)
        ix = it-50:it+50;
        c = c+1;
        xc(it,:) = xcorr(spk(ix)-mean(spk(ix)),'coeff');        
    end;
    
end;

figure;
imagesc(t,-100:100,xc');