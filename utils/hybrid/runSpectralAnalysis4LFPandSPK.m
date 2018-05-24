function [spectralData] = runSpectralAnalysis4LFPandSPK(data_lfp, data_lfp2, data_spk,Fs)

movingwin1 = [0.5 0.01];

TW = 5;
params1                  = [];
params1.tapers           = [TW , 2*TW-1];
params1.pad              = 0;
params1.Fs               = Fs;
params1.err              = 0;
params1.trialave         = 0;
params1.fpass            = [1 30];

TW = 5;
params2                  = [];
params2.tapers           = [TW , 2*TW-1];
params2.pad              = 0;
params2.Fs               = Fs;
params2.err              = 0;
params2.trialave         = 1;
params2.fpass            = [1 30];

%movingwin2 = [0.25 0.01];

% TW = 10*0.25;
% params2                  = [];
% params2.tapers           = [TW , 2*TW-1];
% params2.pad              = 0;
% params2.Fs               = Fs;
% params2.err              = 0;
% params2.trialave         = 0;
% params2.fpass            = [30 200];
tx1 = [];
fx1 = [];
[SgmLFP]   = cell(1,length(data_lfp));
%[SgmLFP2]  = cell(1,length(data_lfp));
for it = 1:length(data_lfp)
    [SgmLFP{it},tx1,fx1] = mtspecgramc(data_lfp2{it},movingwin1,params1);
    %[SgmLFP2{it},tx2,fx2] = mtspecgramc(data_lfp{it},movingwin2,params2);
end;

% [   SgmSPK ]  = cell(1,length(data_spk));
% %[   SgmSPK2]  = cell(1,length(data_spk));
% [      R   ] =  cell(1,length(data_spk));
% %[      R2  ] =  cell(1,length(data_spk));
% for it = 1:length(data_spk)
%     [SgmSPK{it},tx1,fx1,R{it}]   = mtspecgrampb(data_spk{it},movingwin1,params1);
%     %[SgmSPK2{it},~,~,R2{it}] = mtspecgrampb(data_spk{it},movingwin2,params2);
% end;

if (length(data_lfp) == 1) && (length(data_spk) > 1)    
    data_lfp = repmat(data_lfp,[1 length(data_spk)]);
end;

[   Cgrm   ]  = cell(1,length(data_spk));
%[   Cgrm   ]  = cell(length(SgmLFP),length(data_spk));
%[   Cgrm2  ]  = cell(length(SgmLFP),length(data_spk));
%for it = 1:length(SgmLFP)
for it = 1:length(data_spk)
    [Cgrm{it},~, ~,~,~,tx1,fx1] = cohgramcpb(data_lfp{it},data_spk{it},movingwin1,params2);
    %[Cgrm{it,jt}] = cohgramcpb(data_lfp{it},data_spk{it},movingwin1,params1);
    %[Cgrm2{it,jt}] = cohgramcpb(data_lfp{it},data_spk{jt},movingwin2,params2);
end;
%end;

%[   Cgrm3  ]  = cell(length(data_spk),length(data_spk));
%[   Cgrm4  ]  = cell(length(data_spk),length(data_spk));
% for it = 1:length(data_spk)
%     parfor jt = 1:length(data_spk)
%         [Cgrm3{it,jt}] = cohgramcpb(data_spk{it},data_spk{jt},movingwin1,params1);
%         [Cgrm4{it,jt}] = cohgramcpb(data_spk{it},data_spk{jt},movingwin2,params2);
%     end;
% end;

%%
spectralData.SgmLFP = SgmLFP;
%spectralData.SgmLFP2 = SgmLFP2;
%spectralData.SgmSPK = SgmSPK;
%spectralData.SgmSPK2 = SgmSPK2;
spectralData.Cgrm = Cgrm;
%spectralData.Cgrm2 = Cgrm2;
%spectralData.Cgrm3 = Cgrm3;
%spectralData.Cgrm4 = Cgrm4;
spectralData.tx = tx1;
%spectralData.tx2 = tx2;
spectralData.fx = fx1;
%spectralData.fx2 = fx2;
%spectralData.R = R;
%spectralData.R2 = R2;
