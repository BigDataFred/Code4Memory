function [Tlogfile] = set_log4tuning(p2d,pID)
%%

files = dir([p2d,'*.txt']);

Tlogfile = {files(:).name}';