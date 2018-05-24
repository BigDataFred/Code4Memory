function [exp,sesh] = get_EM_rec_info( rpath)

exp = dir(rpath);

chck = regexp({exp.name},'EM');
chck2 = [];
for it = 1:length(chck)
    chck2(it) = ~isempty(chck{it});
end;

exp = {exp(find(chck2)).name};
sesh = cell(1,length(exp));
for st = 1:length(exp)
    
    rpath2 = [rpath,exp{st},filesep];
    
    x = dir(rpath2);
    
    chck = regexp({x(:).name},'\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}');
    chck2 = [];
    for it = 1:length(chck)
        chck2(it) = ~isempty(chck{it});
    end;
    
    [x] = {x(find( chck2 == 1) ).name};% these are the seshions for each patient
    
    c=0;sel_ix=[];
    for it = 1:length(x)
        if regexp(x{it},'\d{4}') ==1
            c =c+1;
            sel_ix(c) = it;
        end;
    end;
    
    sesh{st} = x(sel_ix);
end;