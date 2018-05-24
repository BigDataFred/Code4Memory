function [xc,fr,lag,cnt] = time_resolved_xcorr(dt,lag,trlIx,spkDat1,spkDat2)

trl = cell(1,length(spkDat1));
ts  = cell(1,length(spkDat1));

trl{1} = spkDat1.trial{1};
ts{1} = spkDat1.time{1};
ts{1} = ts{1}.*1e3;

if ~isempty(spkDat2)
    trl{2} = spkDat2.trial{1};
    ts{2} = spkDat2.time{1};
    ts{2} = ts{2}.*1e3;
else
    trl{2} = spkDat1.trial{1};
    ts{2} = spkDat1.time{1};
    ts{2} = ts{2}.*1e3;
end;

%%
selIx = find(dt>=0e3);

%% loop over trials
        
xc = zeros(lag*2+1,length(dt));
cnt = 0;
fr = zeros(length(trlIx),1);
for jt = 1:length(trlIx)
    
    %fprintf([num2str(jt),'/',num2str(length(trlIx))]);
    
    ix1 = find( trl{1} == trlIx(jt));% find spike times corresponding to trial
    ix2 = find( trl{2} == trlIx(jt));% find spike times corresponding to trial
    
    if ( ~isempty(ix1) ) && ( length(ix1) > 1 ) && ( ~isempty(ix2) ) && ( length(ix2) > 1 )
        
        ts2 = ts{1}(ix1);% readout of spike times for trial of interest
        [spk1] = histc(ts2,dt);% do the binning of the spike times
        
        ts2 = ts{2}(ix2);% readout of spike times for trial of interest
        [spk2] = histc(ts2,dt);% do the binning of the spike times
        
        fr(jt) = sum(spk1(selIx))/(length(selIx)/1e3);
                                
        if fr(jt) >5
            
            cnt = cnt+1;    ix2 = [];
                        
            parfor kt = 1:length( dt ) % loop over bins
                
                dum1 = spk1;dum2 = spk2;
                
                if (kt>lag) && (kt < length( dum1 )-lag)
                    ix2 = kt-lag/2:kt+lag/2;% shift window for autocorrelation in 1 ms steps
                    dum = xcorr(dum1(ix2),dum2(ix2));% compute autocorrelation
                    xc(:,kt) = xc(:,kt) +dum';% sum autocorr across trials for each bin
                end;
                
            end;
        end;
        
    end;
    %fprintf('\n');
end;

xc =  xc./cnt;% normalize by the number of trials                
lag = -lag:lag;