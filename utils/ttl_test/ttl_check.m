function [event] = ttl_check(event,basetrig,Fs)
% %%
% addpath('~rouxf/MATLAB/toolboxes/fieldtrip-20160309/');
% ft_defaults;
% addpath('~rouxf/MATLAB/utils/logfile_readers/');    
% %%
% p2d = '/media/rouxf/My Passport/2016-06-03_15-07-24_trigger_test/';
% %%
% [hdr]   = ft_read_header(p2d);
% [event] = ft_read_event(p2d);
%%
x =[];
x(:,1) = [event(:).value];
x(:,2) = [event(:).sample];
%%
ttl =zeros(size(x,1),1);
ttlx =zeros(size(x,1),1);

k = 0;
for it = 2:size(x,1)-1
    
    if ( sign(diff(x(it:it+1))) == -1 ) && ( ( x(it-1) == 0 ) || ( x(it+1) == 128 ) ) && ( x(it,1) ~= 128 )
        
        k = k+1;
        ttl(k) = x(it,1);
        ttlx(k) = it;
        
    end;
    
    if ( sign(diff(x(it:it+1))) == 0 ) && ( x(it+1) == 128 )
        
        k = k+1;
        ttl(k) = x(it,1);
        ttlx(k) = it;
        
    end;
%     
%     if it == 37
%         return;
%     end;
    
    if ( sign(diff(x(it:it+1))) == -1 ) && ( x(it-1) ~= 0 ) && ( x(it-1,1) ~= 128 ) && ( x(it+1) == 0 ) && ( x(it,1) ~= 128 )
        
        %if abs(diff([ttl(k) x(it,1)])) ==1
        if (x(it,2)-x(it-1,2))/Fs >2
        if ( x(it,1) > basetrig )
            k = k+1;
            ttl(k) = x(it,1);
            ttlx(k) = it;
        end;
        end;
        %end;
        
    end;
    
end;
ttl(k+1:end) = [];
ttlx(k+1:end) = [];
% %%
% figure;
% subplot(121);
% plot(x(ttlx,1),'b.');axis tight;
% subplot(122);
% plot(abs(diff(x(ttlx,1))),'r');axis tight;
%%
event = event(ttlx);