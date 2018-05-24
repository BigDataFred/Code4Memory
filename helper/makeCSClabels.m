function [chanID] = makeCSClabels()
%%
c = {};
c{1} = 'CSC-LA';
c{2} = 'CSC-LM';
c{3} = 'CSC-LP';
c{4} = 'CSC-RA';
c{5} = 'CSC-RM';
c{6} = 'CSC-RP';

ch = cell(length(c),8);
for it = 1:8
    for jt = 1:length(c)
        ch(jt,it) = {[c{jt},num2str(it)]};
    end;
end;

chanID = [];
for it = 1:size(ch,1)
    chanID = [chanID ch(it,:)];
    
end;
chanID = chanID';