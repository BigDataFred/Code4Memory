function [refDat] = referenceDat(dat, chanID, chanList,refMethod)

[refDat] = dat;
idx = 0;
for jt = 1:length(refDat.trial)    
    refDat.trial{jt} = zeros(length(chanList)-length(chanID),length(refDat.time{1}));
end;
refDat.label = cell(1,length(chanList)-length(chanID));

for it = 1:length(chanID)
    
    switch refMethod
        case 'commonAvg'
            
            fprintf(['referencing ',chanID{it},' to common average']);
            
            selIx = find(strcmp( chanID(it), chanList));
            
            for jt = 1:length(dat.trial)
                x = dat.trial{jt}(selIx,:);
                commonAVG = ones(size(x,1),1)*mean(x,1);
                refDat.trial{jt}(selIx,:) = dat.trial{jt}(selIx,:) - commonAVG;
            end;
            
        case 'bipolar'
            
            fprintf(['referencing ',chanID{it},' to bipolar montage']);
            
            selIx = find(strcmp( chanID(it), chanList));
            idx = idx(end)+1:idx(end)+length(selIx)-1;
            
            for kt = 1:length(selIx)-1
                for jt = 1:length(dat.trial)
                    x = dat.trial{jt}(selIx(kt:kt+1),:);                    
                    refDat.trial{jt}(idx(kt),:) = x(1,:) - x(2,:);
                end;
                refDat.label{idx(kt)} = [dat.label{selIx(kt)},'-',dat.label{selIx(kt+1)}];
                refDat.label{idx(kt)}(regexp(refDat.label{idx(kt)},' '))=[];
            end;
            
        case 'ica'
            
            fprintf(['decomposing data using ica to find optimal combination']);
            
            cfg                     = [];
            cfg.method              = 'runica';
            
            [refDat] = ft_componentanalysis( cfg , dat );
    end;
    
    fprintf('\n');
end;


return;