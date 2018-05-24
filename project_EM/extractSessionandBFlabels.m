function [sesh,BFlabel,seshLab,dataMode] = extractSessionandBFlabels(pID,expMode,expSel,seshSel)

sesh = [];
BFlabel= [];
seshLab= [];
dataMode= [];

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

if ~isempty(sesh)
    f=0;cnt = 0;
    while f<1
        cnt = cnt+1;
        x = dir([rpath,sesh{cnt},'/lfp_dat/']);% case we're in Bham
        x(1:2) = [];
        if ~isempty(x)
            f=1;
        end;            
    end;
    
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
    
    if isempty(BFlabel)%case we're in Amsterdam
        
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
        if ~isempty(BFlabel)
            dataMode = 'AMS';
        else
            error('unknown dataFormat detected.');
        end;
    else
        dataMode = 'BHX';
    end;
    
    %% extract the data and time of the selected recording
    if ~isempty(seshSel)
        seshLab = sesh{seshSel};
        seshLab(regexp(seshLab,'_')) = [];
    else
        seshLab = cell(1,length(sesh));
        for it = 1:length(sesh)
            seshLab(it) = sesh(it);
            seshLab{it}(regexp(seshLab{it},'_')) = [];
        end;
    end;
end;