function [evdat,evix] = read_tuning_log(p2f,fname)
%%
if nargin == 0
    p2f = '/home/rouxf/Data/EM/Log/P01/SE1/';
    fname = 'LogFile_concept_tuning_P01.txt';
end;
%%
fid = fopen([p2f,fname]);
%%
dat = textscan(fid,'%s');
dat = dat{:};
h = dat(4:11);
%%
dat(1:11) = [];
%%
evdat = reshape(dat,[9 length(dat)/9 ])';
%%
id = unique(evdat(:,4));
evix = cell(1,length(id));
for it = 1:length(id)
    
    [evix{it}] = find(strcmp(evdat(:,4),id{it}));
    
end;