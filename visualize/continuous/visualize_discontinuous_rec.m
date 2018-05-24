%%
f       = double(params.hdr.FirstTimeStamp);

ttlt    = (event.timestamp(:)-(f*g))./1e6/60;
datat   = (params.timestamps(:)-(f*g))./1e6/60';

fst = (f-(f*g))/1e6/60;

%ttlt = ttlt(ttl_idx);

figure;
subplot(141);
hold on;
plot(ttlt,2*ones(1,length(ttlt)),'ko');
plot(datat,ones(1,length(datat)),'r.');
plot(fst,1,'g*');
axis tight;
ylim([-1 4]);

max(ttlt) > max(datat)

dx = diff(datat);
ix = find(dx == max(dx));
sel = find(ttlt >= datat(ix+1));

subplot(142);
plot(dx);
axis tight;
subplot(143);
hold on;
plot(datat,ones(1,length(datat)),'r.');
plot(datat(ix+1),1,'ko');
axis tight;
subplot(144);
f = datat(ix+1).*1e6.*60+f;
ttlt    = (event.timestamp(sel)-(f*g))./1e6/60;
datat   = (params.timestamps(ix+1:end)-(f*g))./1e6/60';

fst = (f-(f*g))/1e6/60;

hold on;
plot(ttlt,2*ones(1,length(ttlt)),'ko');
plot(datat,ones(1,length(datat)),'r.');
plot(fst,1,'g*');
axis tight;
ylim([-1 4]);