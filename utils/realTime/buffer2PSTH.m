%%
connect2Cheetah;
%%
for jt = 1:100
    nbufferRecords = 1;
    bufferSize = .03;
    
    [buffer1] = get_Cheetah_buffer(nbufferRecords,bufferSize,ObjNames{7});

    x1 = buffer1.V{1};
    x1 = x1(:);
    
    ts1 = buffer1.ts{1};
    ts1 =ts1(:);
    ts1 = ts1./1e6;
    
    if any(x1<=-2000)
        break
    end;
end;

nbufferRecords = 1;
bufferSize = 1;

[buffer2] = get_Cheetah_buffer(nbufferRecords,bufferSize,ObjNames(1:4));
x2 = buffer2.V{1};
x2 = x2(:);

ts2 = buffer2.ts{1};
ts2 = ts2(:);
ts2 = ts2./1e6;

figure;
plotyy(ts1,x1,ts2,x2);
%axis tight;
%xlabel('Time (s)');
%ylabel('Amplitude (\muV)');

%%
disconnectCheetah;