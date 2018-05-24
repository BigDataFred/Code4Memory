 %%
ntrlIC = size(ICs1{4},2)/length(tIx) ;
for ot = 8%1:8
for pt = 1:8
%thetaSig = AVG{4};
thetaSig = ICs1{4}(ot,:);
thetaSig = reshape(thetaSig,[length(tIx) ntrlIC]);

gammaSig = ICs2{4}(pt,:);
gammaSig = reshape(gammaSig,[length(tIx) ntrlIC]);
%gammaSig = lfpDat.LFPseg{26};

 pbins2 = 0:20:360;
 phi = zeros(length(trlSel)*length(tIx),length(pfoi));
 amp = zeros(length(trlSel)*length(tIx),length(afoi));
 for jt = 1:ntrlIC
     
     fprintf([num2str(jt),'/',num2str(length(trlSel))]);
     ix = (jt-1)*length(tIx)+1:(jt-1)*length(tIx)+length(tIx);
     
     x = thetaSig(:,jt)';          
     d = zeros(length(tIx),length(pfoi));
     parfor kt = 1:length(pfoi)
         ft = filtfilt(fbp{kt},1,[fliplr(x) fliplr(x) x fliplr(x) fliplr(x)]);
         ft = hilbert(ft);
         %ft = ft(nsmp+1:2*nsmp);
         %ft = ft(tIx);
         ft = ft(2*length(tIx)+1:3*length(tIx));
         d(:,kt) = (angle(ft)+pi).*(180/pi);
     end;
     phi(ix,:) = d;
     
     x = gammaSig(:,jt)';     
     d = zeros(length(tIx),length(afoi));
     for kt = 1:length(afoi)
         ft = filtfilt(fba{kt},1,[fliplr(x) fliplr(x) x fliplr(x) fliplr(x)]);
         ft = hilbert(ft);
         %ft = ft(nsmp+1:2*nsmp);
         %ft = ft(tIx);
         ft = ft(2*length(tIx)+1:3*length(tIx));
         d(:,kt) = abs(ft).^2;
     end;
     amp(ix,:) = d;
     
     fprintf('\n');
 end;
 clear d;
 
 X = zeros(length(afoi),length(pfoi),length(pbins2));
 for kt = 1:length( pfoi )
     p = phi(:,kt);
     for lt = 1:length( afoi )
         a = amp(:,lt);
         parfor mt = 1:length(pbins2)-1
             X(lt,kt,mt) = mean(a(p>=pbins2(mt) & p < pbins2(mt+1)));
         end;
         X(lt,kt,end) = X(lt,kt,1);
         X(lt,kt,:) = X(lt,kt,:)./sum(X(lt,kt,:));
     end;
 end;
 H = -squeeze(sum(X.*log(X),3));clear X;
 n = length(pbins2);
 MI = (log(n)-H)./log(n);clear H;
 
 figure;pcolor(pfoi,afoi,MI);shading interp;colormap jet;%caxis([0 5e-3]);
end;
end;