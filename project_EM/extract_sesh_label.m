function [sesh] = extract_sesh_label(r_path)

[sesh] = dir(r_path);

sel = [];c = 0;
for it = 1:length( sesh )
    
    [chck] = regexp(sesh(it).name,'\w{4}-\w{2}-\w{2}_\w{2}-\w{2}-\w{2}','match');        
    if ~isempty(chck)
        c = c+1;
        sel(c) = it;
    end;
    
end;

[sesh] = sesh(sel);