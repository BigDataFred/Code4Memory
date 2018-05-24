% %% open pool of workers
% if isempty(gcp('nocreate'))
%     parpool(6,'SpmdEnabled',false);
% end;

%%
basepath = '/media/rouxf/rds-share/iEEG_DATA/MICRO/';
pIDs = dir(basepath);
pIDs(1:2) = [];

%%
for kt = 6%1:length(pIDs)
    
    pIDs(kt).name
    
    p2d = [basepath,pIDs(kt).name,filesep];
    chck = dir(p2d);
    chck(1:2) = [];
    for yt = 1:length(chck)
        if ~isempty(regexp(chck(yt).name,'EM'))
            sel = yt;
        end;
    end;
    p2d = [p2d,chck(sel).name,filesep];
    
    fn = dir([p2d,'low_freq_spectrum_',pIDs(kt).name,'.mat']);
    load([p2d,fn.name]);
    
   for jt = 1:length(pow1)
        
        if ~isempty(pow1{jt})
            n = ceil(length(pow1{jt}.label)/8);
            lm = [min(min(pow1{jt}.powspctrm)) min(min(pow2{jt}.powspctrm));...
                max(max(pow1{jt}.powspctrm)) max(max(pow2{jt}.powspctrm))];
            
            lm = [min(lm(:)) max(lm(:))];
            
            a =[];
            figure;
            for it = 1:length(pow1{jt}.label)
                
                subplot(n,8,it);
                a(it) = gca;
                hold on;
                plot(pow1{jt}.freq,pow1{jt}.powspctrm(it,:))
                plot(pow2{jt}.freq,pow2{jt}.powspctrm(it,:),'r');
                axis tight;
                title(pow1{jt}.label(it));
                xlabel('freq');
                ylabel('pow');
            end;
            %set(a,'YLim',lm);
        end;
        
    end;
end;