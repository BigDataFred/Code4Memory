function [sortedSpikes] = segmentSpikeData(sortedSpikes,pre,post,events,lfpTime)    

if nargin == 0
    pre  = 5;
    post = 7;
end;

for jt = 1:length(sortedSpikes)
    
    lfpTime;
    [spkTime] = lfpTime(sortedSpikes{jt}.newSpikeTimes);
    sortedSpikes{jt}.newSpikeTimes = spkTime.*1e6;%convert back from s to micro-s
    sortedSpikes{jt}.SpikeTimesSeg = [];
    sortedSpikes{jt}.trl = [];
    sortedSpikes{jt}.assignedClusterSeg = [];
    sortedSpikes{jt}.oriIx = [];
    
    for it = 1:size(events,1)
        [ix] = find( spkTime >= (events(it,1)-pre) & spkTime <= (events(it,1)+post) );
        sortedSpikes{jt}.oriIx = [sortedSpikes{jt}.oriIx ix];
        sortedSpikes{jt}.SpikeTimesSeg= [sortedSpikes{jt}.SpikeTimesSeg spkTime(ix)-events(it,1)];
        sortedSpikes{jt}.trl = [sortedSpikes{jt}.trl it*ones(1,length(ix))];
        sortedSpikes{jt}.assignedClusterSeg = [sortedSpikes{jt}.assignedClusterSeg sortedSpikes{jt}.assignedCluster(ix)];
    end;
    
end;
    