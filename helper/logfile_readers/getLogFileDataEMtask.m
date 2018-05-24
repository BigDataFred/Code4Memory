function  [logDat]= getLogFileDataEMtask(pID,p2Nlxdata,Nlxdat,sesh)

%% get the logfile data
p = [];
p.pID = pID;
p.sesh = sesh;
p.p2Nlxdata = p2Nlxdata;
p.Nlxdat = Nlxdat;

[sdat] = create_session_dat(p);

params1   = [];
params1.p = sdat.p2lf;%path
params1.fn = sdat.lfn;%filename
params1.ntrl = 49;%number of trials
params1.ncol = 9;% number of columns in logfile

%same as above
params2 = [];
params2.p = sdat.p2lf;
params2.fn = sdat.lfn;
params2.ntrl = 49;
params2.ncol = 12;

%read logfile for Encoding
[LogDat1] = getNewLogDataEM( params1, 'ENC' );
%read logfile for Retrieval
[LogDat2] = getNewLogDataEM( params2, 'RET' );

[RTs] = str2double( LogDat1.log(:,end) );%get RTs
RTs = RTs;

%find hits, misses and stim categories
cat = LogDat1.log(:,3:4);
c = {};
for jt = 1:size( cat,1 )
    
    c{jt} = [ cat{jt,1}(1) cat{jt,2}(1) ];% reads out info from logdata
    
end;

ix    = {};
ix{1} = find( strcmp(c,'fp') );% indices corresponding to face-place tr
ix{2} = find( strcmp(c,'pp') );% indices corresponding to place-place tr
ix{3} = find( strcmp(c,'ff') );% indices corresponding to face-face tr

ix{4} = find( sum( [str2double(LogDat2.log(:,5:6))] ,2) ==2);% both images were correctly remebered
ix{5} = find( sum( [str2double(LogDat2.log(:,5:6))] ,2) ==1);% only 1/2 images were correctly remebered
ix{6} = find( sum( [str2double(LogDat2.log(:,5:6))] ,2) ==0);%no images were correctly remembered

ix_readme = { 'face-place' 'place-place' 'face-face' 'correct-both' 'miss-one' 'miss-both' };% save info for later
LogDat_readme = { 'LogDat1:Enc', 'LogDat2:Ret' };

logDat.LogDat1 = LogDat1;
logDat.LogDat2 =LogDat2;
logDat.ix =ix;
logDat.ix_readme =ix_readme;
logDat.LogDat_readme =LogDat_readme;
logDat.RTs =RTs;
