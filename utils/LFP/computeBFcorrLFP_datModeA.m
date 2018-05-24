function [LFPr] = computeBFcorrLFP_datModeA(dat,corrType,plt)
% extract the median LFP

if isstruct(dat)% case data is formatted for fieldtrip
    f=1;
    nchan = length(dat.label);
    ntrl = length(dat.trial);
else
    f =2;
    nchan = length(dat);
    [nsamp, ntrl] = size(dat{1});
end;
[LFPr] = zeros(ntrl,nchan,nchan);

parfor kt = 1:ntrl
    dat;
    for it = 1:nchan
        if f ==2
            x = dat{it}(:,kt);
        elseif f==1
            x = dat.trial{kt}(it,:)';
        end;
        for jt = 1:nchan
            if f==2
                y = dat{jt}(:,kt);
            elseif f==1
                y = dat.trial{kt}(jt,:)';
            end;
            LFPr(kt,it,jt) = corr(x,y,'Type',corrType);
        end;
    end;
end;

if (nargin>1) && (strcmp(plt,'y'))
    figure;
    imagesc(squeeze(mean(LFPr,1)));axis xy;
    caxis([-1 1]);
    xlabel('Channel #');
    ylabel('Channel #');
end;
