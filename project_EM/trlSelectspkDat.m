function [spkDat] = trlSelectspkDat(spkDat,trlENC)

for jt = 1:length( spkDat )
    for kt = 1:length( spkDat{jt}.unit )
        
        if ~isempty( spkDat{jt}.unit{kt} )
            
            selIx = find(ismember( spkDat{jt}.trial{kt},trlENC )==1);
            
            spkDat{jt}.timestamp{kt}            = spkDat{jt}.timestamp{kt}(selIx);
            spkDat{jt}.waveform{kt}             = spkDat{jt}.waveform{kt}(:,:,selIx);
            spkDat{jt}.unit{kt}                 = spkDat{jt}.unit{kt}(selIx);
            spkDat{jt}.time{kt}                 = spkDat{jt}.time{kt}(selIx);
            spkDat{jt}.trial{kt}                = spkDat{jt}.trial{kt}(selIx);
        end;
        
    end;
end;