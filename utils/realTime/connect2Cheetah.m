%%
restoredefaultpath;
addpath('C:\\toolbox\MatlabImportExport_v6.0.0\');     
addpath('C:\\Netcom\Matlab_M-files\');
%addpath(genpath('C:\\Users\rouxf-admin\Desktop\Stream Channel\'));    
%%
serverName = '10.2.2.30';
succeeded = NlxConnectToServer(serverName)
%%
if NlxAreWeConnected()
    NlxSetApplicationName('Matlab real time app');
    [succeeded, ObjNames, ObjTypes] = NlxGetDASObjectsAndTypes();
else
    NlxDisconnectFromServer();
    return;    
end;
%%
succeeded = zeros(length(ObjNames),1);
for it = 1:length(ObjNames)
    [succeeded(it)] = NlxOpenStream(ObjNames{it});
end;
