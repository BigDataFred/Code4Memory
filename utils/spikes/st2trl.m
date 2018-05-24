function [trlst] = st2trl(st,trlt,Fs,dt)

dt = abs(dt);

if length(dt)==1
    dt = [dt dt];
end;

trlst = cell(length(trlt),1);
for it = 1:length(trlt)
    
    dx = linspace(trlt(it)-dt(1),trlt(it)+dt(2),Fs);
        
    trlst{it} = st(st >= min(dx) & st <= max(dx));
        
    dx = dx - trlt(it);
    trlst{it} = trlst{it} - trlt(it);
    
end;

return;