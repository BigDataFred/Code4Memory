function [delIx,pct,chLabIED] = IEDdetector(LFPdat,chanLab,toi)

% cfg                     = [];
% cfg.lpfilter            = 'yes';
% cfg.lpfiltord           = 4;
% cfg.lpfreq              = 30;
% cfg.lpfilttype          = 'but';
% 
% [LFPdat] = ft_preprocessing( cfg , LFPdat );

cfg                     = [];
cfg.channel             = chanLab;
cfg.keeptrials          = 'yes';

[LFPdat] = ft_timelockanalysis( cfg, LFPdat);
chLabIED = LFPdat.label;

cfg                     = [];
cfg.latency             = [toi(1) toi(2)];

[LFPdat] = ft_selectdata( cfg , LFPdat );

ERPimg = LFPdat.trial;

[delIx,pct] = IEDdetector2(ERPimg);


% 
% trl = {};
% %figure;
% for jt = 1:size(LFPdat.trial,2)
%     dum = [];
%     for it = 1:size(LFPdat.trial,1)
%         
%         
%         
%         x = squeeze(LFPdat.trial(it,jt,:));
%         y = x.^2;
%         z = (y-mean(y))./std(y);
% %                 subplot(211);
% %                 plot(x);
% %                 subplot(212);
% %                 hold on;
% %                 plot(z);
% %                 plot([1 length(z)],[median(z)+15*iqr(z) median(z)+15*iqr(z)],'r--');
% %                 pause;
% %                 clf;
%         thr = z> (median(z)+20*iqr(z));
%         if any(thr)
%             dum=  [dum it];
%         end;
%         
%     end;
%     trl{jt} = dum;
% end;
% 
% pct = [];
% for it = 1:length(trl)
%     pct(it) = length([trl{it}])/size(LFPdat.trial,1);
% end;