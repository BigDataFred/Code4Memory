function [datac] = glmcleaning(comp,ic_idx,varargin)

datac = comp;

%figure;
c2 = 0;
for zt = 1:length(varargin)
    zt
    %%
    rc = zeros(size(comp.sampleinfo,1),2);
    n = length(comp.trial{1});
    for it = 1:length(comp.trial)
        rc(it,:) = [(it-1)*n+1 (it-1)*n+n];
    end;
    
    chck = [];
    I = zeros(size(varargin{1}{zt}.sampleinfo,1),1);
    for it = 1:size(varargin{1}{zt}.sampleinfo,1)
        
        m = [varargin{1}{zt}.sampleinfo(it,1) varargin{1}{zt}.sampleinfo(it,2)];
        
        chck(:,1) = sign(rc(:,1)-m(1));
        chck(:,2) = sign(rc(:,2)-m(2));
        
        x =  [find(sum(chck,2)==0)];% find(sum(chck,2)==-2)];
        if isempty(x)
            x =  max(find(sum(chck,2)==-2));
        end;
        
        if isempty(x)
            I(it) = 0;
        else
            I(it) = x;
        end;
        
    end;
    %%
    id = unique(I);
    id(id==0) = [];
    c = 0;
    n = length(comp.trial{1});
    X = zeros(4,length(comp.trial{1})*length(comp.trial));
    X(1,:) = 1;%constant term
    for it = 1:length(id)
        
        ix =  find(I == id(it));
        
        m = [varargin{1}{zt}.sampleinfo(ix,:)];
        
        for jt = 1:size(m,1)
            c = c+1;
            ix2 = m(jt,1):m(jt,2);
            %cf = eval([num2str(log10(max(abs(varargin{1}{zt}.avg))))]);
            cf = round(abs(log10(max(abs(varargin{1}{zt}.avg)))));
            %if sign(abs(cf)-1)==-1 && sign(cf)==-1;cf = 10^abs(cf);else;cf =1;end;                        
            
            X(2,ix2) = c*ones(1,length(ix2));%artefact number 
            
            X(3,ix2) = it*ones(1,length(ix2));%trial number
            
            X(4,ix2) = varargin{1}{zt}.avg.*10^cf;% average waveshape
            
            %X(4,ix2) = varargin{1}{zt}.trl(c,:);% single trial waveshape
            %X(4,ix2) = X(4,ix2).*10^cf;
            
            %X(5,ix2) = mean(varargin{1}{zt}.trl(c,:))+std(varargin{1}{zt}.trl(c,:),0,2)*randn(1,length(ix2));% random noise
            %X(5,ix2) = X(5,ix2).*10^cf;
            
        end;
    end;
    del_idx = find(X(2,:)==0);
    sel_idx = setdiff(1:size(X,2),del_idx);
    
    % M  = mean(X(2,:));
    % SD = std(X(2,:),0,2);
    % for jt = 1:size(X,1)
    %     for it = 1:length(del_idx)
    %         X(jt,del_idx(it)) = M+SD*randn(1,1);
    %     end;
    % end;
    
    %ic_idx(zt)_idx = 1:size(X,2);
    %%
    Y = zeros(length(comp.label),length(comp.trial{1})*length(comp.trial));
    
    idx = 1:length(comp.trial{1});
    for it = 1:length(comp.trial)
        Y(:,idx) = comp.trial{it};
        idx = idx + length(comp.trial{it});
    end;
    
    %%
    %[B] = glmfit(X([5],sel_idx)',Y(sel_idx)');
    for it = 1:length(comp.label)
    [B] = regress(Y(it,sel_idx)',X(:,sel_idx)');
    %%
    yp = zeros(1,size(Y,2));
    
    for lt = 1:length(B)
        yp(sel_idx) = B(lt)*X(lt,sel_idx);
    end;
    
    %yp(sel_idx) = glmval(B,X([5]),'identity');%
    
    %%
    dc = Y(it,:)-yp;
            
    idx = 1:length(comp.trial{1});
    n = length(idx);
    for jt = 1:length(comp.trial)
        datac.trial{jt}(it,:) = dc(idx);
        
        idx = idx+n;
    end;
%     %%
%     cfg = [];
%     cfg.method = 'mtmfft';
%     cfg.pad = 'maxperlen';
%     cfg.output = 'pow';
%     cfg.taper = 'dpss';
%     cfg.tapsmofrq = .5;
%     cfg.channel = comp.label(ic_idx(zt));
%     
%     [pow1] = ft_freqanalysis(cfg,comp);
%     
%     [pow2] = ft_freqanalysis(cfg,datac);
    
    %%  
%     s = length(ic_idx);
%     
%     c2 = c2+1;
%     subplot(s,3,c2);
%     plot(Y);
%     hold on;
%     plot(yp,'r');
%     %plot(yp,'c');
%     
%     c2 = c2+1;
%     subplot(s,3,c2);
%     hold on;
%     plot(Y);
%     plot(Y-yp,'k');
%      
%     c2 = c2+1;
%     subplot(s,3,c2);
%     hold on;
%     plot(pow1.freq,10*log10(pow1.powspctrm));
%     plot(pow2.freq,10*log10(pow2.powspctrm),'r');
%     xlim([0 100]);
    end;
end;
%%



return;