%         %% clean line noise
%         s = 15;
%         [step] = round(length(dataSamples)/(s*Fs));
%         
%         ix = 1:Fs*s;
%         [dum] = dataSamples;
%         nf = [];
%         for jt = 1:step
%             
%             fprintf([num2str(jt),'/',num2str(step)]);
%             
%             flag = 0;
%             %while flag ==0
%             
%             T = length(dum(ix))/Fs;
%             W = 0.1;%1/T;
%             TW = round(T*W);
%             k = round(2*TW-1);
%             params                  = [];
%             params.Fs               = Fs;
%             params.pad              = 8;
%             params.fpass            = [40 60];
%             params.trialave         = 0;
%             params.tapers           = [TW k];
%             
%             [S,f] = mtspectrumc( dum(ix)', params );
%             S = 20*log10(S);
%             S = (S-mean(S))./std(S);
%             
%             figure;plot(f,S);hold on;
%             
%             S2 = S.*(S>2);%plot(f,S2,'r');
%             dx = sign(gradient(S2));
%             ix1 = find(dx==1);
%             ix2 = find(dx==-1);
%             ix1 = ix1(find(gradient(ix1)>1));
%             ix2 = ix2(find(gradient(ix2)>1));
%             
%             if length(ix1) ~= length(ix2)
%                 if ix2(1)<ix1(1);ix2(1) = [];end;
%                 chck = [];
%                 for yt = 1:length(ix1);
%                     for zt= yt:length(ix2);
%                         chck(yt,zt) = ix1(yt)>ix2(zt);
%                     end;
%                 end;
%                 [ix1,ix2]= find(chck);
%                 ix2(ix2) = [];
%                 
%             end;
%             
%             if ~isempty(ix1) && ~isempty(ix2)
%                 d = pdist(f([ix1 ix2]),'euclidean');
%                 z = linkage(d,'ward');
%                 c = cluster(z,'maxclust',3);
%                 
%                 %                 figure;
%                 %                 hold on;
%                 %                 plot(f(ix1(find(c==1))),f(ix2(find(c==1))),'bs');
%                 %                 plot(f(ix1(find(c==2))),f(ix2(find(c==2))),'rs');
%                 %                 plot(f(ix1(find(c==3))),f(ix2(find(c==3))),'gs');
%                 
%                 pIx = [1 2;1 3;2 3];
%                 d= [];
%                 for zt = 1:size(pIx,1)
%                     d(zt) = abs(diff(round([mean(f(ix1(find(c==pIx(zt,1))))) mean(f(ix1(find(c==pIx(zt,2)))))])));
%                 end;
%                 d = d>=1;
%                 
%                 ix4 = pIx(find(d==0),:);
%                 
%                 ix3 = pIx(find(d==1),:);
%                 ix3 = ix3(:);
%                 ix3( ismember(ix3,c(find( ismember(c, ix4)))) ) = [];
%                 ix3 = unique(ix3);
%                 
%                 f2 = [];
%                 if ~isempty(ix4)
%                     x = [ix1( find( ismember(c, ix4))) ix2( find( ismember(c, ix4)))];
%                     n = [];
%                     for zt = 1:size(x,1)
%                         n(zt) = length(x(zt,1):x(zt,2));
%                     end;
%                     [~,sIx] = sort(n);
%                     cIx = length(sIx)-1:length(sIx);
%                     cIx = cIx(cIx>0);
%                     f2 = mean(mean(f(x(sIx(cIx),:))));
%                 end;
%                 
%                 if ~isempty(ix3)
%                     for zt = 1:length(ix3)
%                         x = [ix1( find( ismember(c, ix3(zt)))) ix2( find( ismember(c, ix3(zt))))];
%                         n = [];
%                         for yt = 1:size(x,1)
%                             n(yt) = length(x(yt,1):x(yt,2));
%                         end;
%                         [~,sIx] = sort(n);
%                         cIx = length(sIx)-1:length(sIx);
%                         cIx = cIx(cIx>0);
%                         f2 = [f2;mean(mean(f(x(sIx(cIx),:))))];
%                     end;
%                 end;
%                 f2 = sort(f2);
%                 f2 = f2(f2>47 & f2 <53);
%                 
%                 if ~isempty(f2)
%                     %                     x1 = [];
%                     %                     x2 = [];
%                     %                     t = ix./Fs;
%                     %                     for it = 1:length(f2)
%                     %                         x1(:,it) = sin(2*pi*f2(it).*t);
%                     %                         x2(:,it) = cos(2*pi*f2(it).*t);
%                     %                     end;
%                     %
%                     %                     X = [x1 x2 sum(x1,2) sum(x2,2) x1(:,1).*x1(:,2) x2(:,1).*x2(:,2) x1(:,1).*x2(:,1) x1(:,2)+x2(:,2) x1(:,2).*x2(:,2) x1(:,1)+x2(:,2) x1(:,1).*x2(:,2) x1(:,2)+x2(:,1) x1(:,2).*x2(:,1)];
%                     %
%                     %                     Y =dum(ix)';
%                     %                     Y = Y-mean(Y);
%                     %
%                     %                     b = regress(Y,[ones(size(X,1),1) X]);
%                     %
%                     %                     yp = [ones(size(X,1),1) X]*b;
%                     %
%                     %                     %
%                     %                     %                 figure;
%                     %                     %                 subplot(121);
%                     %                     %                 plot(t,Y)
%                     %                     %                 xlim([0 .1])
%                     %                     %                 hold on
%                     %                     %                 plot(t,yp,'r');
%                     %                     %                 subplot(122);
%                     %                     %                 hold on;
%                     %                     %                 plot(t,Y)
%                     %                     %                 plot(t,Y-yp,'r');xlim([0 .1])
%                     %                     %
%                     %                     dum(ix) = dum(ix) -yp';
%                     %
%                     
%                     for zt = 1:length(f2)
%                         [dum(ix)] = CleanLineNoise(dum(ix),'Fs',Fs,'noiseFreq',f2(zt),'windowSize',s);
%                     end;
%                     [S,f] = mtspectrumc( dum(ix)', params );
%                     S = 20*log10(S);
%                     S = (S-mean(S))./std(S);
%                     plot(f,S,'r');
%                 end;
%             end;
%             
%             %                 ix1 =[find(f<=48) find(f >=51)];
%             %                 [~,m] = max(S(ix1));
%             %                 ix1 = ix1(m);
%             %                 ref1 = S(ix1);
%             %
%             %                 ix2 = find(f>=48 & f <=51);
%             %                 [~,m] = max(S(ix2));
%             %                 ix2 = ix2(m);
%             %                 ref2 = S(ix2);
%             %
%             %                 if ref2 > ref1
%             %                     nf = [nf f(ix2)];
%             %                     [dum(ix)] = CleanLineNoise(dum(ix),'Fs',Fs,'noiseFreq',f(ix2),'windowSize',s);
%             %                 else
%             %                     flag = 1;
%             %                 end;
%             %
%             %end;
%             ix = ix+s*Fs;
%             if max(ix) > length(dataSamples)
%                 ix = ix(1):length(dataSamples);
%                 s = length(ix)/Fs;
%             end;
%             fprintf('\n');
%         end;
%         
%         [dataSamples] = dum;