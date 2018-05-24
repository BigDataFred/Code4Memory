function [buffer] = get_Cheetah_buffer(nbufferRecords,bufferSize,objName)

V = cell(nbufferRecords,1);
itimeStampArray = cell(nbufferRecords,1);
k=0;
for it = 1:nbufferRecords+1
    
    pause(bufferSize);
    
    parfor jt = 1:length(objName)
        [succeeded, dataArray, timeStampArray, channelNumberArray, samplingFreqArray, numValidSamplesArray, numRecordsReturned, numRecordsDropped ] = NlxGetNewCSCData(objName{jt});
    end;
    
    if it >=2
        k = k+1;
        [V{k}] = convertAD2V(dataArray,timeStampArray,numRecordsReturned);
        
        [timeStampArray,V{k}] = sort_NLx_records(timeStampArray,V{k});        
        [itimeStampArray{k}] = interpTimeStampArray(double(timeStampArray),double(unique(samplingFreqArray)));        
    end;
end;

buffer.V = V;
buffer.ts = itimeStampArray;