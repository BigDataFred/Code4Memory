function [LogDat] = getNewLogDataTune(f,p)
%%
if nargin == 0
    f.p2logf = '/home/rouxf/Data/EM/Log/P01/';
    f.logf = 'P01_24-May-2015_00_00_0000000_log_tune.txt';
    p.ncols = 8;
end;
%%
fid = fopen([f.p2logf,f.logf],'r');

dat = textscan(fid,'%s');
[LogDat.hdr] = dat{1}(1:3);

dat{1}(1:3) = [];

idx = 1:8;
f = 0;
x = {};
while f<1
    y = dat{1}(idx)';      
    x = [x;y];
    idx = idx+8;
    if max(idx)>length(dat{1})
        f=1;
    end;
end;

idx=1:8;f = 0;x = [];while f<1;x = [x;dat{1}(idx)'];idx = idx+8;if max(idx)>length(dat{1});f=1;end;end;
LogDat.dat =x;
[LogDat.dat] = reshape(dat{1},[p.ncols length(dat{1})/p.ncols])';

if ~isequal(x,LogDat.dat )
    error('LogFile mismatch detected');
end;

LogDat.dat(1,:) = [];

[LogDat.RT] = str2double(LogDat.dat(:,end));

[LogDat.stimID] = str2double(LogDat.dat(:,1));

[LogDat.nS,LogDat.xS] = hist(LogDat.stimID,[min(LogDat.stimID):max(LogDat.stimID)]);
% if max(LogDat.nS) ~= min(LogDat.nS)
%     error('frequency of events is out of range');
% end;

LogDat.ID = unique(LogDat.stimID);
LogDat.n = size(LogDat.dat,1);