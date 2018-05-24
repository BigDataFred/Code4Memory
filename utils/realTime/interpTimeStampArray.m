function [itimeStampArray] = interpTimeStampArray(timeStampArray,Fs)

itimeStampArray = cell(length(timeStampArray),1);
for jt = 1:length(timeStampArray)-1
    itimeStampArray{jt} = linspace(timeStampArray(jt),timeStampArray(jt+1)-1/Fs,512);
end;
itimeStampArray{jt+1} = linspace(timeStampArray(jt+1),timeStampArray(jt+1)+(1/Fs)*512,512);

itimeStampArray = [itimeStampArray{:}];

return;
