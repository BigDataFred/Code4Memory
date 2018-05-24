function [trigSampEnc,trigSampRet,recT] = reconstructTTLtimingFromLog(ixON,ixOFF,trigDat,LogDat1,LogDat2)
%%
retIx = find(diff([ixON' ixOFF'],[],2)./trigDat.fsample>1.9);

te = trigDat.time{1}(ixON(retIx(1)));
retT = str2double(LogDat2.log(1,10));

encT = str2double(LogDat1.log(:,5));
dumT = zeros(length(encT),1);

c = 0;
for it = length(encT):-1:1
    c = c+1;
    
    dumT(c) =retT - encT(it);
    
end;

recT = te-dumT;

trigSamp = zeros(1,length(dumT));
x = 1:length(trigDat.trial{1});
for it = 1:length(recT)
    trigSamp(it) = nearest(x,recT(it)*trigDat.fsample);
end;
trigSamp = fliplr(trigSamp);

trigSampEnc = trigSamp;
trigSampRet = ixON(retIx);

figure;
plot(trigDat.time{1},trigDat.trial{1},'b-');
hold on;
plot(recT,zeros(1,length(dumT)),'go');
plot(trigDat.time{1}(trigSampEnc),trigDat.trial{1}(trigSampEnc),'r.');
plot(trigDat.time{1}(trigSampRet),trigDat.trial{1}(trigSampRet),'k.');

return;