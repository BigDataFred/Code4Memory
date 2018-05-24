function [spectralDat] = spectralAnalysisEM2(lfpDat,chanIx,trlIx,pre,tIx,lF,hF,trialave)

%%
addpath(genpath('/home/rouxf/tbx/chronux_2_11/'));
addpath(genpath('/home/rouxf/tbx/eeglab14_1_1b/functions/sigprocfunc/'));
addpath(genpath('/home/rouxf/prj/Bham/code/mcode/utils/'));

%% recruit workers for parallel computing
if isempty(gcp('nocreate'))
    parpool(36,'SpmdEnabled', false)
end;

%%
[paramsTF1,paramsTF2,movingwin1,movingwin2,paramsSP1,paramsSP2] = setParameters4spectralAnalysis(lfpDat.Fs,fix(length(tIx)/lfpDat.Fs),trialave);

%% compute power spectra
spectralDat = cell(length(chanIx),1);
parfor it = 1:length( chanIx )
    fprintf([num2str(it),'/',num2str(length( chanIx ))]);
    spectralDat{it} = struct;
        
    if isempty(trlIx)
        trlIx2 = 1:size(lfpDat.LFPseg{chanIx(it)},2);
    else
        trlIx2 = trlIx{it};
    end;
        
    if (lF == true)  && ( ~isempty(trlIx2) )
        if ~isempty(tIx)
            [S1,f1,Serr1] = mtspectrumc( gradient(lfpDat.LFPseg{chanIx(it)}(tIx,trlIx2)')', paramsSP1);
        else
            S1 = [];f1 = [];Serr1 = [];
        end;
        [S3,t3,f3] = mtspecgramc( gradient(lfpDat.LFPseg{chanIx(it)}(:,trlIx2)')', movingwin1, paramsTF1);
        t3 = t3-pre;
        spectralDat{it}.fx1 = f1;
        spectralDat{it}.Sx1 = S1;
        spectralDat{it}.Sxe = Serr1;
        spectralDat{it}.txx1 = t3;
        spectralDat{it}.fxx1 = f3;
        spectralDat{it}.Sxx1 = S3;
    end;
    if (hF == true) && ( ~isempty(trlIx2) )
        if ~isempty(tIx)
            [S2,f2,Serr2] = mtspectrumc( gradient(lfpDat.LFPseg{chanIx(it)}(tIx,trlIx2)')', paramsSP2);
        else
            S2 = [];f2 = [];Serr2 = [];
        end;
        [S4,t4,f4] = mtspecgramc( gradient(lfpDat.LFPseg{chanIx(it)}(:,trlIx2)')', movingwin2, paramsTF2);
        t4 = t4-pre;
        spectralDat{it}.fx2 = f2;
        spectralDat{it}.Sx2 = S2;
        spectralDat{it}.Sxe = Serr2;
        spectralDat{it}.txx2 = t4;
        spectralDat{it}.fxx2 = f4;
        spectralDat{it}.Sxx2 = S4;
    end;
    fprintf('\n');
end;
	  
            
