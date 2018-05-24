function [dat1,ixON,ixOFF] = runIEDdetection(dat1,dat2,dt,plt)

%%
ixON  = cell( length( dat1.trial ) , length(dat2.label) );
ixOFF = cell( length( dat1.trial ) , length(dat2.label) );

%%
if strcmp(plt,'y');figure;end;

%%
for it = 1:length( dat1.trial )
    fprintf([num2str(it),'/',num2str(length( dat1.trial ))]);
    %%
    [dum1] = dat1;
    for jt = 1:length( dum1.label)
        dum1.label{jt}(end+1) = '1';
    end;
    
    [dum2] = dat2;
    for jt = 1:length( dum2.label)
        dum2.label{jt}(end+1) = '2';
    end;
    %%
    [avg] = ft_timelockanalysis([],dum2);
    
    scf = 70/median(mean(avg.avg,2));
    
    %%
    thr = 3.5;
    Fs = 1/diff(fliplr(dum2.time{1}(1:2)));
    
    for jt = 1:length(dum2.label)
        [ixON{it,jt},ixOFF{it,jt}] = detectIED(dum1.trial{it}(jt,:),dum2.trial{it}(jt,:).*scf,thr,dt,plt);
        if strcmp(plt,'y')
            pause;
            clf;
        end;
    end;
    
    fprintf('\n');
%     %%
%     [dum] = ft_appenddata( cfg , dum1 , dum2);
%     
%     cfg                     = [];
%     cfg.keeptrials          = 'yes';
%     
%     [dum] = ft_timelockanalysis( cfg , dum );
% 
%     dum = rmfield(dum,'var');
%     dum.dimord = 'chan_time';
%     dum.avg = squeeze(dum.trial(it,:,:));
%     dum = rmfield(dum,'trial');
%          
%     a = zeros( 1 , length( dum.label ));
%     m = zeros( 2 , length( dum.label ));
%     for jt = 1:length( dum.label )
%         subplot(length(dum.label),1,jt);
%         a(jt) = gca;
%         m(:,jt) = [min(dum.avg(jt,:)) max(dum.avg(jt,:))];
%         plot(dum.time,dum.avg(jt,:),'k');
% %         if ~isempty(ixON{jt})
% %             subplot(length(dum.label),1,jt);
% %             hold on;
% %             plot(dum.time(ixON{jt}),dum.avg(jt,ixON{jt}),'r^','MarkerFaceColor','r');
% %             plot(dum.time(ixOFF{jt}),dum.avg(jt,ixOFF{jt}),'g^','MarkerFaceColor','g');
% %             axis off;
% %         end;
%     end;
%     axis(a,'tight');
%     axis(a,'off');
%     set(a,'YLim',[min(min(m)) max(max(m))]);     
%     for jt = 1:length(a)
%         text(a(jt),dum.time(10),max(max(m)),dum.label);
%     end;
%     pause;
%     clf;
    
end;
