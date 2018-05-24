%%
succeeded = zeros(length(ObjNames),1);
for it = 1:length(ObjNames)
    [succeeded(it)] = NlxCloseStream(ObjNames{it});
end;
%%
%Disconnect from the server.
NlxDisconnectFromServer();
fprintf('Disconnected from %s.\r', serverName);