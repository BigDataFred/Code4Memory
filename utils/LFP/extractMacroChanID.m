
function [chanID] = extractMacroChanID(chanLabels)

ix = regexp(chanLabels,'\d{1}');

chanID = cell( length( ix ) ,1 );
for it = 1:length( ix )
    
    chanID(it) = {chanLabels{it}(1:ix{it}-1)};       
    chanID{it}(regexp(chanID{it},' ')) = [];
    
end;

return