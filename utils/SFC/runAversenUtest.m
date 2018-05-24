function [dz,vdz,Adz,fBinz] = runAversenUtest(data1c1,data1c2,data2c1,data2c2,params,pval,mode)
%%
[N,~] = size(data1c1);

nfft = 2^(nextpow2(N)+params.pad);
df = 1/(nfft/params.Fs);
fBinz = params.fpass(1):df:params.fpass(2);

[J1c1] = preData4AversenTestc(data1c1,params);
[J1c2] = preData4AversenTestc(data1c2,params);


if mode ==1
    [dz,vdz,Adz]=two_group_test_spectrum_corrected(J1c1,J1c2,pval);
else
    if ~isempty(data2c1) && ~isempty(data2c2)
        [J2c1] = preData4AversenTestpb(data2c1,params);
        [J2c2] = preData4AversenTestpb(data2c2,params);
        [dz,vdz,Adz]=two_group_test_coherence_corrected(J1c1,J2c1,J1c2,J2c2,pval);
    end;
end;

    
