function [LFPdataClean] = cleanLFPfromLineNoise(LFPdata,Fs,step,s)

ix = 1:Fs*s;
[dum] = LFPdata;

for jt = 1:step % loop over the number of subsequent time windows
    
    fprintf([num2str(jt),'/',num2str(step)]);
    
    flag = 0;
    while flag<1
        
        T = length(dum(ix))/Fs;
        W = 1/T;
        TW = round(T*W);
        k = round(2*TW-1);
        params                  = [];
        params.Fs               = Fs;
        params.pad              = 8;
        params.fpass            = [40 60];
        params.trialave         = 0;
        params.tapers           = [TW k];
        
        % estimate noise frequency
        [S,f] = mtspectrumc( dum(ix)', params );
        S = 20*log10(S);
        S = (S-mean(S))./std(S);
        
        %ref = mean([S(f<=45);S(f>55)]);
        selIx = find(f>=47 & f <=51);
        [m,mIx] = max(S(selIx));
        
        if ( m > iqr([S(f<=45);S(f>55)])*2 )
            % remove line noise through subtraction of fitted line
            [dum(ix),~] = CleanLineNoise(dum(ix),'Fs',Fs,'noiseFreq',f(selIx(mIx)),'windowSize',s);
        else
            flag =1;
            [S2] = mtspectrumc( dum(ix)', params );
            S2 = 20*log10(S2);
            S2 = (S2-mean(S2))./std(S2);
        end;
    end;
    
    ix = ix+s*Fs;
    if max(ix) > length(dum)
        ix = ix(1):length(dum);
        s = length(ix)/Fs;
    end;
    fprintf('\n');
    
end;

[LFPdataClean] = dum; 