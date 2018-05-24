function [trlst] = st2event(st,et,Fs,dt)

dt = abs(dt);

if length(dt)==1
    dt = [dt dt];
end;

trlst = cell(length(et),1);
for it = 1:length(et)
    
    d = linspace(et(it)-dt,et(it)+dt,Fs);
        
    trlst{it} = st(st >= min(d) & st <= max(d));
        
    d = d - et(it);
    trlst{it} = trlst{it} - et(it);
    
end;

return;