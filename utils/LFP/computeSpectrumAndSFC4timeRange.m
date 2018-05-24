function [spectrumLFP,spectrumSFC] = computeSpectrumAndSFC4timeRange(data_lfp , data_spk, timeIdx,trlIdx, Fs)

%% LFP-power and SFC specrta

if isempty(timeIdx)
    timeIdx = 1:size(data_lfp{1},1);
end;
if isempty(trlIdx)
    trlIdx = 1:size(data_lfp{1},2);
end;

spectrumLFP = cell(1,length(data_lfp));
%spectrumSFC = cell(length(data_lfp),length(data_spk));
spectrumSFC = cell(1,length(data_lfp));

%%
for it = 1:length(data_lfp)
   
    sel_data1 = data_lfp{it}(timeIdx,trlIdx);
    [spectrumLFP{it}] = computeLFPspectrum(sel_data1,Fs);
    
    if ~isempty(data_spk)
        %for jt = 1:length(data_spk)
            %sel_data2 = data_spk{jt}(timeIdx,trlIdx);
            %[spectrumSFC{it,jt}] = computeSFCspectrum(sel_data1, sel_data2, Fs);
            sel_data2 = data_spk{it}(timeIdx,trlIdx);
            [spectrumSFC{it}] = computeSFCspectrum(sel_data1, sel_data2, Fs);
        %end;
    end;
    
end;
clear sel_data*;