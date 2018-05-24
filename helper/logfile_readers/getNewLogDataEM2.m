function [dE,dR,b] = getNewLogDataEM2(params)

params.ncols1 = 9;
params.ncols2 = 12;

%%
fid = fopen([params.p,params.fn],'r');

[dat] = textscan(fid,'%s');

%%
sel = [];
chck = regexp(dat{1},'ENC');
for it = 1:length(chck)
    if ~isempty(chck{it})
        sel = [sel it];
    end;
end;

sel2 = [];
chck = regexp(dat{1}(sel),'\d{1}');
for it = 1:length(chck)
    if ~isempty(chck{it})
        sel2 = [sel2 it];
    end;
end;

[selE] = sel(sel2);

%%
sel = [];
chck = regexp(dat{1},'RET');

for it = 1:length(chck)
    if ~isempty(chck{it})
        sel = [sel it];
    end;
end;

sel2 = [];
chck = regexp(dat{1}(sel),'\d{1}');
for it = 1:length(chck)
    if ~isempty(chck{it})
        sel2 = [sel2 it];
    end;
end;

[selR] = sel(sel2);

%%
dE = [];
dR = [];
b = [];
for it = 1:length(selE)
    
    x = dat{1}(selE(it)+1:selR(it)-1)';
    x = reshape(x,[params.ncols1 length(x)/params.ncols1])';
    dE = [dE;x];
    
    b = [b size(x,1)];
    
    if it < length(selE)
        x = dat{1}(selR(it)+1:selE(it+1)-5)';
    else
        x = dat{1}(selR(it)+1:length(dat{1})-4)';
    end;
    x = reshape(x,[params.ncols2 length(x)/params.ncols2])';
    dR = [dR;x];
end;
    

