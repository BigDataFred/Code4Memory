function [Nlxdat] = set_Nlx4tuning(p2d,pID)
%%

files = dir([p2d,'*']);
files([1:2 find(strcmp({files(:).name}','res'))]) = [];

Nlxdat = {files(:).name}';
