function sdat = create_session_dat(params)
%%
%sdat.p2Nlxf = [filesep,'media',filesep,'rouxf',filesep,'rds-share',filesep,'Archive',filesep,'MICRO',filesep,piD,filesep,'fvsPEM',filesep];
%x = dir(sdat.p2Nlxf);
%sdat.p2Nlxf = [sdat.p2Nlxf,x(sesh+2).name,filesep];
%sdat.p2Nlxf = params.p2Nlxdata;

ix2 = regexp(params.p2Nlxdata,'EM')-1;
ix1 = max(regexp(params.p2Nlxdata(1:ix2),'/'))+1;
expName = params.p2Nlxdata(ix1:ix2);

sdat.p2lf  = [filesep,'media',filesep,'rouxf',filesep,'rds-share',filesep,'Archive',filesep,'MICRO',filesep,params.pID,filesep,'log',filesep,'EM',filesep];

sdat.lfn = dir(sdat.p2lf);
pattern = [params.pID,'_',expName];
pattern2 = [params.pID(2:end),'_',expName];

chck1 = [];
chck2 = [];
for it = 1:length(sdat.lfn)
    chck1(it) = ~isempty(regexpi(sdat.lfn(it).name,pattern));
    chck2(it) = ~isempty(regexpi(sdat.lfn(it).name,pattern2));
end;
chck = chck1+chck2;
sdat.lfn= sdat.lfn(chck~=0);

x = {};f1 = 0; f2 = 0;
for it = 1:length(sdat.lfn);
    chck1 =  regexp(sdat.lfn(it).name,'\d{2}-\w{3}-20\d{2}','match');
    chck2 =  regexp(sdat.lfn(it).name,'20\d{2}-\w{3}-\d{2}','match');
    
    if ~isempty(chck1)
        f1 =1;
        x(it)=chck1;
    else
        f2 = 1;
        x(it)=chck2;
    end;
    
end;

if f2 ~=0 % account for inconsistent format of log files
    x2 = {};
    for it = 1:length(x)
        x2{it} = [x{it}(max(regexp(x{it},'-'))+1:end) x{it}(min(regexp(x{it},'-')):max(regexp(x{it},'-'))) x{it}(1:min(regexp(x{it},'-'))-1)];
    end;
    
    dum = datenum(x2);
    
    [~,ix] = sort(dum);
    x2 = x2(ix);
    
    x ={};
    for zt = 1:length(x2)
        dum = datevec(x2{zt});
        dum = dum(1:3);
        x{zt} = [];
        for it = 1:length(dum)
            if dum(it)<10
                x{zt} = [x{zt} '0',num2str(dum(it)),'-'];
            else
                x{zt} = [x{zt} num2str(dum(it)),'-'];
            end;
            
        end;
        x{zt}(end) = [];
    end;
end;

sesh_sel = datenum(datestr(params.Nlxdat{params.sesh}(1:regexp(params.Nlxdat{params.sesh},'_')-1)));

sesh_all = datenum(x);
%[~,ix] = sort(sesh_all);
[ix] = find(sesh_all == sesh_sel);

% if sesh_all(params.sesh) ~= sesh_sel
%     error('logfile assignment does not match');
% end;

% if ~strcmp(params.Nlxdat{params.sesh}(1:regexp(params.Nlxdat{params.sesh},'_')-1),x{params.sesh});
%     error('logfile assignment does not match');
% end;

sdat.lfn = sdat.lfn(ix(ismember(ix,params.sesh))).name;

%sdat.lfn = sdat.lfn(params.sesh).name;

return;