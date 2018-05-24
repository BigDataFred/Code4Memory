function visualizeUnits4manualSelection(wvfStats,instFR,spkRaster,dt,saveParams)
%%
h = figure;set(h,'Visible','off');
if ~isempty(spkRaster)
    subplot(4,5,[4 5 9 10]);
    hold on;
    plot(linspace(0,2,64),wvfStats{1}.wvf,'r');
    plot(linspace(0,2,64),mean(wvfStats{1}.wvf,2),'k');
    axis tight;
    ylabel('Amplitude [\muV]');
    xlabel('Time [ms]');
    
    % subplot(4,5,[1 2 6 7 11 12]);
    % hold on;
    % for jt = 1:size( spkRaster ,2 )
    %     r = spkRaster{1,jt};
    %     t = jt*ones(1,length(r));
    %     r = [r;r];
    %     t = [t-.5;t+.5];
    %     line(r,t,'Color','k');
    %     r = [];
    % end;
    % plot([0 0],[0 size(spkRaster,2)+1],'Color','r');
    % plot([2e3 2e3],[0 size(spkRaster,2)+1],'Color','r');
    % xlim([-.5 5].*1e3);
    % ylim([0 size(spkRaster,2)+1]);
    % ylabel('Trial #');
    % xlabel('Time [ms]');
    
    subplot(4,5,[14 15]);
    bar(0:250,wvfStats{1}.nIsi);
    axis tight;
    ylabel('Count');
    xlabel('ISI [ms]');
    
    subplot(4,5,[19 20]);
    plot(wvfStats{1}.lag,wvfStats{1}.xc);
    axis tight;
    ylabel('Coincidences');
    xlabel('Lag [ms]');
    
    subplot(4,5,[16 17]);
    y = instFR{1};
    hold on;
    plot(dt(1:end-1),y(1:end-1),'k-s','LineWidth',3);
    plot([0 0],[0 max(y)+1],'Color','r');
    plot([2e3 2e3],[0 max(y)+1],'Color','r');
    axis tight;xlim([-.5 5].*1e3);
    ylabel('Spikes/s');
    xlabel('Time [ms]');
end;
saveas(gcf,[saveParams.savepath,saveParams.savename],'fig');
close;

