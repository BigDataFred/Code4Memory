function [LFPseg] = segmentLFPdata(LFPsig,trl,nsmp)

LFPseg = cell( 1 , length(LFPsig) );
for jt = 1:length(LFPsig)
    
    fprintf([num2str(jt),'/',num2str(length(LFPsig))]);
    
    x = zeros(size(trl,1),nsmp);
    for it = 1:size(trl,1)
        x(it,:) = LFPsig{jt}(trl(it,1):trl(it,2));
    end;
    LFPseg{jt} =  x';
    
    fprintf('\n');
end;
clear trl;