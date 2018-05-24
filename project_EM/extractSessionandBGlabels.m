function [sesh,BFlabel,seshLab] = extractSessionandBGlabels(pID,expMode,expSel,seshSel)

[rpath] = ['/media/rouxf/rds-share/iEEG_DATA/MICRO/',pID,'/',expMode{expSel},'/'];

x = dir(rpath);
c = 0;chck = [];
for it = 1:length(x)
    if ~isempty(regexp(x(it).name,'\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}'))
        if x(it).isdir
            c = c+1;    chck(c) = it;
        end;
    end;
end;
[sesh] = {x(chck).name}';

x = dir([rpath,sesh{1},'/lfp_dat/']);
pattern1 = 'CSC_';
pattern2 = '_\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}';
c = 0;  BFlabel = {};
for it = 1:length(x)
    ix = [regexp(x(it).name,pattern1)+4:regexp(x(it).name,pattern2)-1];
    if ~isempty(ix)
        c = c+1;    BFlabel(c) = {x(it).name(ix)};  BFlabel{c}(regexp(BFlabel{c},'\d{1}')) = [];
    end;
end;
[BFlabel] = unique(BFlabel);

if ~isempty(BFlabel)
    
    pattern1 = 'lfp_data_';
    pattern2 = '_\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}';
    c = 0;  BFlabel = {};
    for it = 1:length(x)
        ix = [regexp(x(it).name,pattern1)+9:regexp(x(it).name,pattern2)-1];
        if ~isempty(ix)
            c = c+1;    BFlabel(c) = {x(it).name(ix)};  BFlabel{c}(regexp(BFlabel{c},'\d{1}')) = [];
        end;
    end;
    [BFlabel] = unique(BFlabel);
    
end;

%% extract the data and time of the selected recording
seshLab = sesh{seshSel};
seshLab(regexp(seshLab,'_')) = [];