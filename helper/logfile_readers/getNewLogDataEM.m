function [LogDat] = getNewLogDataEM(params,mode)

%%
fid = fopen([params.p,params.fn],'r');

dat = textscan(fid,'%s');

%%
k = 0;
chck1 = regexp(dat{:},[mode,'\d{1}']);

idx1 = [];
for it =1:length(chck1)
    
    if ~isempty(chck1{it})
                
        k = k+1;
        idx1(k,1) = it;        
        
    end;
    
end;

if ~isempty(idx1)
    idx1 = idx1+1;
end;

%%
k = 0;
switch mode
    case 'ENC'
        chck2 = regexp(dat{:},'RET\d{1}');
    case 'RET'
        chck2 = regexp(dat{:},'BreakBeg\d{1}');
        
end;

idx2 = [];
for it =1:length(chck2)
    
    if ~isempty(chck2{it})
                
        k = k+1;
        idx2(k,1) = it;        
        
    end;
    
end;

if ~isempty(idx2)
    idx2 = idx2-1;
end;

%%
if isempty(idx2)
    idx2 = length(dat{:});
end;

%%
n1 = length(idx1);
n2 = length(idx2);
if n1 > n2
    %idx1 = idx1(1:n2); % change added 29/08/2017 by FRO
    idx2(end+1) = length(dat{:});
elseif n2 > n1
    idx2 = idx2(1:n1);
end;

%%
d = cell(params.ntrl,params.ncol);
k = 0;
idx3 = zeros(length(idx1),1);
n = zeros(length(idx1),1);
for it = 1:length(idx1)
    
    x = dat{1}(idx1(it):idx2(it));
    idx = 1:params.ncol;
    
    for jt = 1:length(x)/params.ncol            
        k = k+1;
        if any(strcmp(x(idx),'NaN'))
            ix = find(strcmp(x(idx),'NaN'));
            idx(end-1:end) = [];
            
            d(k,:) = [x([idx(1:ix(2))])' 'NaN' 'NaN' x(idx(ix(2)+1:end))'];
        elseif sum(str2double(x(idx)) ==0) ==2
            if sum(isnan(str2double(x(idx(end-1:end))))) ~=0
                ix = find(str2double(x(idx))==0);
                idx(end-1:end) = [];                
                d(k,:) = [x([idx(1:ix(2))])' '0' '0' x(idx(ix(2)+1:end))'];
            else
                d(k,:) = x(idx);
            end;
        else
            d(k,:) = x(idx);
        end;
        idx = idx(end)+1:idx(end)+params.ncol;
    end;
    idx3(it) = k;
    n(it) = length(x)/params.ncol;
    
    
end;
d(k+1:end,:) = [];

%%
%if size(d,1) ~= sum(n);error('trial numbers do not match up');end;
%%
LogDat.log = d;
LogDat.idx = zeros(length(idx3),2);
for it =1:length(idx3)
    if it == 1
        LogDat.idx(it,:) = [1 idx3(it)];
    else
        LogDat.idx(it,:) = [idx3(it-1)+1 idx3(it)];
    end;
end;

%%
return;