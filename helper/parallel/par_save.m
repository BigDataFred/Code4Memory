function par_save(savename,v,varargin)

save_data = cell(1,length(varargin));
for it = 1:length(varargin)
    save_data{it} = varargin{it};
end;

if isempty(v)
    eval(['save(''',savename,''',''save_data'')']);
else
    eval(['save(''',savename,''',''save_data'',''',v,''')']);
end;
