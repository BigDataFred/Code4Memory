function [LogDat] = make_log_data_tune(bf)
if nargin == 0

    bf.p2logf = '/home/rouxf/Data/EM/Log/P01/';
    bf.logf = 'P01_25-May-2016_10_13_3612100_log_tune.txt';
    bf.paramf = 'P01_tuning_param_25-May-2016_10:13:44.mat';
    bf.savepath = [bf.p2logf,'LogMat/'];
    
end;
%% get the logfile data
p.ncols = 8;
[LogDat] = getNewLogDataTune(bf,p);
%%
LogDat.cond = LogDat.dat(2:end,3);
LogDat.stimID = str2double(LogDat.dat(2:end,1));
%%
load([bf.p2logf,bf.paramf]);
LogDat.params = params;
%%
sn = [bf.logf(1:regexp(bf.logf,'.txt')-1),'.mat'];
save([bf.savepath,sn],'LogDat');

return;