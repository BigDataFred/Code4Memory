function [tsi] = timeStampinter(ts)

x = 1:512:(length(ts))*512;
x(end+1) = x(end)+512-1;

y = ts;
dt = y(2)-y(1);
y(end+1) = y(end)+dt;

xi = 1:x(end);
yi = interp1(x,y,xi,'linear');


tsi = yi;



if any(isnan(yi))
    error('time stamp data must be integer');
end;