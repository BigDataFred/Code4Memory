function [fr] = compute_firing_rate(pSTH,dt)
    

[fr] = sum(pSTH,1)/(size(pSTH,1)*(dt));

return;