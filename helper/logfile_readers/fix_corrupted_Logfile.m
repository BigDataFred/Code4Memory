function [LogDat] = fix_corrupted_Logfile(fname,p,nreps,npic)
%%
if nargin ==0
    fname.p2logf = [''];
    fname.logf = 'P05_TS01_TS01_log_ctune_26112016_14_48_30.txt';
    nreps = 6;
    npic = 59;
end;

fid = fopen([fname.p2logf,fname.logf]);
data=textscan(fid,'%s');
data = data{1};

LogDat.hdr = data(1:3);

data(1:3) = [];
data(1:p) = [];

nev = npic*nreps;

idx = 1:p;
dat = {};
for it = 1:nev
    try
        x = {};
        x = data(idx);
        
        c=3;
        chck = length(data{idx(c)})==1;
        while ~chck
            c =c+1;
            chck = length(data{idx(c)})==1;
        end;
        f = p-(2+length(idx(c:end)));
        
        if c>3
            idx2 = [idx(1:2) idx(c:end) idx(end)+1:idx(end)+f];
            x = data(idx2);
            x(2) = {[data{idx(2)} data{idx(3)}]};
            idx = idx2(end)+1:idx2(end)+p;
        else
            idx = idx + p;
        end;
        
        dat(it,:) = x';
    catch
    end;
end;

% if ~( nev == size(dat,1) )
%     error(' number of events does not match ');
% end;

LogDat.dat = dat;

[LogDat.RT] = str2double(LogDat.dat(:,end));

[LogDat.stimID] = str2double(LogDat.dat(:,1));

[LogDat.nS,LogDat.xS] = hist(LogDat.stimID,[min(LogDat.stimID):max(LogDat.stimID)]);
% if max(LogDat.nS) ~= min(LogDat.nS)
%     error('frequency of events is out of range');
% end;

LogDat.ID = unique(LogDat.stimID);
LogDat.n = size(LogDat.dat,1);