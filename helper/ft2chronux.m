function [data_lfp,data_spk, data_spk2] = ft2chronux(spkDat,MWdat,BFinf,trlIdx,toi)

%% prep data for Chronux - use for LFP power analysis
if ~isempty(spkDat) 
    [data_spk]          = cell(1,length(spkDat));
    [data_spk2]         = cell(1,length(spkDat)); %format for STA
    cnt = 0;sel = [];
    for jt = 1:length(spkDat)
        [~,~,data_spk{jt},spk2] = prepSPKandLFPdat4Chronux(spkDat,[],[],jt,trlIdx,[toi(1) toi(2)],MWdat.fsample);
        data_spk2{jt} = spk2;
        if ~isempty(data_spk{jt})
            cnt = cnt+1;
            sel(cnt) = jt;
        end;
    end;
    [data_spk] = data_spk(sel);
    [data_spk2] = data_spk2(sel);
end;

if ~isempty(MWdat) 
    [data_lfp]          = cell(1,length(MWdat));
    cnt = 0;sel = [];
    for jt = 1:length(MWdat.label)
        [data_lfp{jt}] = prepSPKandLFPdat4Chronux([],MWdat,jt,[],trlIdx,[toi(1) toi(2)],MWdat.fsample);
        if ~isempty(data_lfp{jt})
            cnt = cnt+1;
            sel(cnt) = jt;
        end;
    end;
    
    [data_lfp] = data_lfp(sel);
end;
data_lfp = data_lfp([BFinf{2,:}]);

return;

%% this section can be used for spike-LFP interpolation, but better to do it
%%before downsampling of data

if ~isempty(spkDat) && ~isempty(MWdat)
    [data_lfp_interp]   = cell(1,length(spkDat));
    cnt = 0; ix = [];
    for it = 1:length(spkDat)
        
        if any(strcmp(spkDat{it}.hdr.label,MWdat.label))
            cnt = cnt+1;
            ix(cnt) = find(strcmp(spkDat{it}.hdr.label,MWdat.label));
            if ~strcmp(spkDat{it}.hdr.label , [BFinf{1,it}]);
                error('index assignment is out of range');
            end;
            [~,data_lfp_interp{cnt}] = prepSPKandLFPdat4Chronux(spkDat,MWdat, BFinf{2,it},it,trlIdx,[toi(1) toi(2)],MWdat.fsample);
            
        end;
        
    end;
    cnt = 0;
    for it = 1:length(ix)
        cnt = cnt+1;
        data_lfp{ix(it)} = data_lfp_interp{cnt};
    end;
    clear data_lfp_interp;
end;

