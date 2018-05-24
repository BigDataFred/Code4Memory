function [LogDat] = getNewLogData(f,p)
%%
if nargin == 0
    f.path2file = '/media/rouxf/My Passport/LOG/';
    f.filename = 'P01_25-May-2016_9_27_5275100_LogFile_EMtask.txt';
    p.ntrl = 100;
    p.ncol = 11;
end;

%%
fid = fopen([f.path2file,f.filename],'r');

dat = textscan(fid,'%s');
%%
k = 0;
chck1 = regexp(dat{:},'RET\d{1}');

idx1 = [];
for it =1:length(chck1)
    
    if ~isempty(chck1{it})
                
        k = k+1;
        idx1(k,1) = it;        
        
    end;
    
end;
%%
k = 0;
chck2 = regexp(dat{:},'BreakBeg\d{1}');

idx2 = [];
for it =1:length(chck2)
    
    if ~isempty(chck2{it})
                
        k = k+1;
        idx2(k,1) = it;        
        
    end;
    
end;
%%
d = cell(p.ntrl,p.ncol);
k = 0;
idx3 = zeros(length(idx1),1);
for it = 1:length(idx1)
    
    x = dat{1}(idx1(it)+1:idx2(it)-1);
        idx = 1:11;
        
        for jt = 1:length(x)/p.ncol
            k = k+1;
            d(k,:) = x(idx);
            idx = idx+11;
        end;
        idx3(it) = k;
end;
%%
LogDat.log = d;
LogDat.idx = zeros(length(idx3),2);
for it =1:length(idx3)
    if it == 1
        LogDat.idx(it,:) = [1 idx3(it)-1];
    else
        LogDat.idx(it,:) = [idx3(it-1) idx3(it)-1];
    end;
end;
%%
return;