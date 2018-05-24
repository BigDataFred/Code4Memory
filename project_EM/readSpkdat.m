function [spkDat,uIx,x,chLab] = readSpkdat(p2SpkDat,spkFiles,trlENC)

chLab = cell(1,length(spkFiles));
spkDat = cell(1,length(spkFiles));
for it = 1:length(spkFiles)
    fprintf([num2str(it),'/',num2str(length(spkFiles))]);
    dat = load([p2SpkDat,spkFiles(it).name]);
    spkDat{it} = dat.save_data{1}{1}{1};
    chLab(it) = {spkDat{it}.hdr.label};
    fprintf('\n');
end;
%extract the spk data that corresponds to the encoding trials
for it = 1:length( spkDat )
    fprintf([num2str(it),'/',num2str(length(spkFiles))]);
    
    if ~isempty( spkDat{it}.unit{1})
        for jt = 1:length( spkDat{it}.unit)
            
            % find indexes that link spike times to trials
            selIx = find( ismember( spkDat{it}.trial{jt}, trlENC ) );%only keep the selected trials
            % keep only selected trials
            spkDat{it}.unit{jt}         = spkDat{it}.unit{jt}(selIx);
            spkDat{it}.trial{jt}        = spkDat{it}.trial{jt}(selIx);
            spkDat{it}.timestamp{jt}    = spkDat{it}.timestamp{jt}(selIx);
            spkDat{it}.time{jt}         = spkDat{it}.time{jt}(selIx);
            spkDat{it}.waveform{jt}     = spkDat{it}.waveform{jt}(:,:,selIx);
        end;
    end;
    
    fprintf('\n');
end;

%% assign an index to each MW and to each cluster on each MW
x = {}; c = 0;  uIx = [];

for it = 1:length( spkDat )% loop over the MW
    
    fprintf([num2str(it),'/',num2str( length( spkDat ) ) ]);
    if ~isempty(spkDat{it}.unit{1})
        for jt = 1:length( spkDat{it}.unit) % loop over the units 

            trl = spkDat{it}.trial{jt};% extract the trial indexes
            st  = spkDat{it}.time{jt}.*1e3;% convert spike times from s to ms
            
            c = c+1;
            uIx(c,:) = [it jt]; % 1st index: MW, second index: unit/cluster
                        
            for kt = 1:length( trlENC )% loop over the encoding trials                
                x{c,kt} = st(trl == trlENC(kt));% use this to make raster plot               
            end;
            clear trl st;
            
        end;
    end;
    fprintf('\n');
    
end;