function [XC,shXC] = computeClusterXcorr(spkDat,toi,lag,trlENC,plt,chLab)
%%
if nargin ==1
    error('you must specify a time range over which to compute the correlations');
end;

dt = toi(1):toi(2);
XC = {}; shXC = {};
p = []; c =0; 
for it = 1:length(spkDat)        
    
    
    trl1 = spkDat{it}.trial{1};
    ts1 = spkDat{it}.time{1}.*1e3;
    delIx = [find(ts1 <toi(1)) find(ts1 > toi(2))];
    trl1(delIx) = [];
    ts1(delIx) = [];    
    
    cnt = 0;
    for jt = 1:length(spkDat)
        if it~=jt
            cnt = cnt+1;
            c = c+1;
            p(c,:) = [it jt];
            trl2 = spkDat{jt}.trial{1};
            ts2  = spkDat{jt}.time{1}.*1e3;
            delIx = [find(ts2 <toi(1)) find(ts2 > toi(2))];
            trl2(delIx) = [];
            ts2(delIx) = [];
                        
            xc = zeros(length(trlENC),length(-lag:lag));
           
            parfor kt = 1:length(trlENC)
                x1 = ts1(trl1 == trlENC(kt));
                x2 = ts2(trl2 == trlENC(kt));
                
                [n1,~] = hist(x1,dt);
                [n2,~] = hist(x2,dt);
                
                xc(kt,:) = xcorr(n1,n2,lag);
            end;
            
            xc2 = zeros(1e3,length(trlENC),length(-lag:lag));
            
            %FIXME
%             parfor st = 1:2e3          
%                 trlENC; ts1; ts2; dt;
%                 sh1 = trlENC(randperm(length(trlENC)));
%                 sh2 = trlENC(randperm(length(trlENC)));
%                 dumxc = zeros(length(trlENC),length(-lag:lag));
%                 for kt = 1:length(sh1)
%                     x1 = ts1(trl1 == sh1(kt));
%                     x2 = ts2(trl2 == sh2(kt));
%                     x1(x1<dt(1)) = [];
%                     x1(x1>dt(end)) = [];
%                     x2(x2<dt(1)) = [];
%                     x2(x2>dt(end)) = [];
%                     [n1,~] = hist(x1,dt);
%                     [n2,~] = hist(x2,dt);
%                     
%                     dumxc(kt,:) = xcorr(n1,n2,lag);
%                 end;    
%                 xc2(st,:,:) = dumxc;
%             end;
            
            XC{it,cnt} = xc;
            shXC{it,cnt} = xc2;
        end;
        
    end;
end;

%%
if strcmp(plt,'y') || strcmp(plt,'yes')
    figure;
    n1 = size(XC,1);
    n2 = size(XC,2);
    cnt = 0;
    for it = 1:size(XC,1)
        for jt = 1:size(XC,2)
            cnt = cnt+1;
            subplot(n1,n2,cnt);
            hold on;
            bar(-lag:lag,mean(XC{it,jt},1));
            %plot(-lag:lag,max(squeeze(mean(sum(shXC{it,jt},1),2)))*ones(1,length(-lag:lag)),'r--');
            axis tight;
            if ~isempty(chLab)
                ylabel(chLab(p(cnt,1)));
                title(chLab(p(cnt,2)));
            end;
        end;
    end;
end;
