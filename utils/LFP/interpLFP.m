function [lfp_interp] = interpLFP(lfp,spk,Fs,method,dt)

if nargin <5
    dt = [0.002 0.006];
end;
dt = dt.*Fs;

N = length(lfp);
lfp_interp = lfp;

if [min(unique(spk)) max(unique(spk))] == [0 1]
    sel = find(spk);
else
    sel = spk;
end;

% figure;
for kt = 1:length(sel)
    
    spk_ix = [];
    spk_ix = sel(kt)-dt(1):sel(kt)+dt(2);
    spk_ix = spk_ix(spk_ix-2 >=1 & spk_ix+2 <=N);
               
    if ~isempty(spk_ix)
        sel2 = [spk_ix(1)-2:spk_ix(1)-1 spk_ix(end)+1:spk_ix(end)+2];
        
        y = lfp(sel2);
        x = sel2;
        
        switch method
            case 'linear'
                yi = interp1(x,y,spk_ix,'linear');
            case 'spline'
                yi = spline(x,y,spk_ix);
        end;
        
        if any(isnan(yi))
            error('NaN detected');
        end;
        
        lfp_interp(spk_ix) = yi;
        
%         subplot(5,5,kt);
%         figure;
%         hold on;
%         plot(spk_ix(1)-100:spk_ix(end)+100,lfp(spk_ix(1)-100:spk_ix(end)+100),'k');
%         plot(spk_ix,yi,'c');
%         plot(spk_ix,yi2,'m');
%         axis tight;box off;
%         if kt == 1
%             legend('raw LFP','linear','spline');
%         end;
%         xlabel('Samples');
%         ylabel('[\muV]');
        
    end;
   
end;

