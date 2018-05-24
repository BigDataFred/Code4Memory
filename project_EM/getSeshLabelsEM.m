function [sesh] = getSeshLabelsEM(rpath)
sesh = dir(rpath);

if ~isempty(sesh)

sesh(1:2) = [];

[chck] = regexp({sesh(:).name},'\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}');

sel = [];
for it = 1:length(chck)
    if ~isempty(chck{it})
        sel = [sel it];
    end;
end;
sesh = sesh(sel);
clear sel chck;

dum = cell(length(sesh),1);
for it = 1:length(sesh);
    dum(it) = {[sesh(it).name]};
end;
sesh = dum;
else
sesh = [];
end;
